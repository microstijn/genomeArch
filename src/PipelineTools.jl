module PipelineTools

using CSV
using DataFrames
using Statistics
using GLM
using StatsModels
using Optim
using CairoMakie

export merge_and_impute_ogt
export optimize_models
export mechanistic_engine
export plot_evolutionary_trajectories

# Phase 0: OGT Database Integration & Taxonomic Imputation
"""
    merge_and_impute_ogt(genarch_file::String, tempura_file::String, output_file::String)

Merges genArch output with the TEMPURA OGT database via NCBI TaxId, and imputes
missing OGT values using a taxonomic roll-up (median of the closest parent clade).
"""
function merge_and_impute_ogt(genarch_file::String, tempura_file::String, output_file::String)
    println("Reading genArch and TEMPURA data...")
    # Read files
    df_genarch = CSV.File(genarch_file) |> DataFrame
    df_tempura = CSV.File(tempura_file) |> DataFrame

    # Standardize taxonomy_id column name for merging
    if "taxonomy_id" in names(df_tempura)
        rename!(df_tempura, :taxonomy_id => :taxId)
    end

    # We will need the taxonomy columns from genarch later for imputation.
    # TEMPURA also has taxonomic info, but we only want to keep its data columns
    # to avoid conflict/duplication with genarch's taxonomy, so we drop tempura's
    # taxonomic columns except taxId.
    cols_to_drop = filter(x -> match(r"^(superkingdom|phylum|class|order|family|genus|species|strain)$"i, x) !== nothing, names(df_tempura))
    df_tempura_clean = select(df_tempura, Not(cols_to_drop))

    # Perform Left Join
    println("Merging datasets via TaxId...")
    df_merged = leftjoin(df_genarch, df_tempura_clean, on=:taxId, makeunique=true)

    println("Imputing missing OGT values via taxonomic roll-up...")
    if !("Topt_ave" in names(df_merged))
        @warn "Column 'Topt_ave' not found in TEMPURA database. Cannot impute OGT."
        CSV.write(output_file, df_merged)
        return df_merged
    end

    # Build reference dictionaries for fast median lookups instead of filtering repeatedly
    function build_clade_medians(df, level)
        if !(string(level) in names(df))
            return Dict{String, Float64}()
        end
        # Filter rows with valid OGT and non-missing level
        valid_rows = filter(row -> !ismissing(row[level]) && !ismissing(row.Topt_ave) && !isnan(row.Topt_ave), df)

        medians_dict = Dict{String, Float64}()
        if nrow(valid_rows) > 0
            grouped = groupby(valid_rows, level)
            for g in grouped
                key_val = string(first(g)[level])
                # Filter out "NA" or empty strings if they crept in
                if key_val != "NA" && key_val != ""
                    medians_dict[key_val] = median(g.Topt_ave)
                end
            end
        end
        return medians_dict
    end

    levels = [:genus, :family, :order, :class, :phylum, :superkingdom]
    clade_medians = Dict(lvl => build_clade_medians(df_merged, lvl) for lvl in levels)

    # Global median fallback
    valid_global = filter(row -> !ismissing(row.Topt_ave) && !isnan(row.Topt_ave), df_merged)
    global_median = nrow(valid_global) > 0 ? median(valid_global.Topt_ave) : NaN

    imputed_topt = Float64[]
    for row in eachrow(df_merged)
        if ismissing(row.Topt_ave) || isnan(row.Topt_ave)
            val = missing
            # Traverse lineages upwards
            for level in levels
                if string(level) in names(df_merged) && !ismissing(row[level])
                    key_val = string(row[level])
                    if haskey(clade_medians[level], key_val)
                        val = clade_medians[level][key_val]
                        break
                    end
                end
            end

            if ismissing(val)
                val = global_median
            end
            push!(imputed_topt, val)
        else
            push!(imputed_topt, Float64(row.Topt_ave))
        end
    end

    df_merged.OGT = imputed_topt

    println("Writing merged and imputed data to $output_file...")
    CSV.write(output_file, df_merged)

    return df_merged
end

# Phase 1: The "Environmentally-Aware" AICc Optimizer
"""
    optimize_models(df::DataFrame)

Calculates Convergent Overlap Density (overlaps per 1000 genes) and compares
three competing models using AICc to determine the best predictor.
"""
function optimize_models(df::DataFrame)
    println("Preparing data for model optimization...")
    # Calculate target variable: Convergent Overlaps per 1000 genes
    # C_overlap_nr represents the count of convergent overlaps.
    if !("C_overlap_nr" in names(df))
        @error "Data missing 'C_overlap_nr' column. Cannot calculate C_overlap_density."
        return nothing
    end

    # Needs total genes. Calculate if not present.
    if !("total_genes" in names(df))
        if "p_gene_nr" in names(df) && "n_gene_nr" in names(df)
            df.total_genes = df.p_gene_nr .+ df.n_gene_nr
        else
            @error "Data missing gene count columns ('p_gene_nr', 'n_gene_nr' or 'total_genes')."
            return nothing
        end
    end

    # Filter out rows with 0 genes just in case
    df_clean = filter(row -> row.total_genes > 0 && !ismissing(row.OGT) && !isnan(row.OGT), df)

    # Genome Size in Mb
    if !("genome_size_mb" in names(df_clean))
        if "genome_size" in names(df_clean)
            df_clean.genome_size_mb = df_clean.genome_size ./ 1_000_000
        elseif "contig_size" in names(df_clean)
            df_clean.genome_size_mb = df_clean.contig_size ./ 1_000_000
        else
            @error "Data missing 'genome_size' or 'contig_size'."
            return nothing
        end
    end

    # Target variable
    df_clean.C_overlap_density = (df_clean.C_overlap_nr ./ df_clean.total_genes) .* 1000

    # Derive is_free_living boolean.
    # In AnalysisTools.jl, host_associated is marked by certain regex.
    # Here we'll recreate a boolean for free-living.
    if "is_free_living" in names(df_clean)
        # Already present
    elseif "environment" in names(df_clean)
        # Re-derive from environment string. Free-living is NOT host-associated.
        df_clean.is_free_living = [ismissing(env) || !occursin(r"host|human|rumen|feces|gut|symbiont|pathogen"i, env) ? 1.0 : 0.0 for env in df_clean.environment]
    elseif "is_host_associated" in names(df_clean)
        df_clean.is_free_living = 1.0 .- df_clean.is_host_associated
    else
        @warn "No environment information found. Assuming all are free-living (1.0) for the sake of Model 3 testing."
        df_clean.is_free_living = ones(nrow(df_clean))
    end

    println("Fitting Model 1 (The Rigid Baseline)...")
    # k = 2 (Intercept + genome_size_mb)
    mod1 = lm(@formula(C_overlap_density ~ genome_size_mb), df_clean)

    println("Fitting Model 2 (The Thermal Hypothesis)...")
    # k = 3 (Intercept + genome_size_mb + OGT)
    mod2 = lm(@formula(C_overlap_density ~ genome_size_mb + OGT), df_clean)

    println("Fitting Model 3 (The Decoupled Lifestyle)...")
    # k = 4 (Intercept + genome_size_mb + OGT + genome_size_mb * is_free_living)
    # The prompt says: "overlaps are driven by Size and OGT, but the rate of overlap accumulation is decoupled by a Lifestyle boolean".
    # This implies an interaction between Size and the Lifestyle boolean.
    mod3 = lm(@formula(C_overlap_density ~ genome_size_mb + OGT + genome_size_mb & is_free_living), df_clean)

    # Function to calculate AICc
    function calc_aicc(model, k, n)
        # Extract log likelihood
        ll = loglikelihood(model)
        aic = -2 * ll + 2 * k
        if n - k - 1 > 0
            aicc = aic + (2 * k * (k + 1)) / (n - k - 1)
        else
            aicc = aic # fallback if sample size is extremely small
        end
        return aicc
    end

    n_obs = nrow(df_clean)

    aicc1 = calc_aicc(mod1, 2, n_obs)
    aicc2 = calc_aicc(mod2, 3, n_obs)
    aicc3 = calc_aicc(mod3, 4, n_obs)

    println("\n--- AICc Model Comparison ---")
    println("Model 1 (Baseline):           AICc = $aicc1")
    println("Model 2 (Thermal Hypothesis): AICc = $aicc2")
    println("Model 3 (Decoupled Lifestyle):AICc = $aicc3")

    models = Dict(1 => (mod1, aicc1, "Model 1 (The Rigid Baseline)"),
                  2 => (mod2, aicc2, "Model 2 (The Thermal Hypothesis)"),
                  3 => (mod3, aicc3, "Model 3 (The Decoupled Lifestyle)"))

    best_mod_idx = argmin([aicc1, aicc2, aicc3])
    best_model, best_aicc, best_name = models[best_mod_idx]

    println("-> The winning model is: $best_name")

    return best_model, df_clean
end

# Phase 2: The Decoupled Mechanistic Engine
"""
    mechanistic_engine(df::DataFrame)

Calculates the predicted C_overlap_density by identifying the critical
intergenic spacing threshold below which a "toxicity penalty" (overlap surge) activates.
"""
function mechanistic_engine(df::DataFrame)
    println("Initializing the Decoupled Mechanistic Engine...")

    if !("C_overlap_density" in names(df))
        if "C_overlap_nr" in names(df) && "total_genes" in names(df)
            df.C_overlap_density = (df.C_overlap_nr ./ df.total_genes) .* 1000
        else
            @error "Missing 'C_overlap_density' and cannot calculate it."
            return nothing
        end
    end

    if !("mean_gap_size" in names(df))
        @error "Missing 'mean_gap_size' (intergenic spacing). Cannot find threshold."
        return nothing
    end

    # We want to model:
    # C_overlap_density = base_rate + penalty(mean_gap_size)
    # penalty(x) = max(0, a * (threshold - x))
    # where base_rate could be a small constant or a function of genome size.
    # To isolate the threshold effect, we optimize parameters (base_rate, a, threshold)
    # to minimize the Sum of Squared Errors (SSE) or a robust loss function.

    # Let's clean data
    valid_data = filter(row -> !ismissing(row.mean_gap_size) && !isnan(row.mean_gap_size) &&
                               !ismissing(row.C_overlap_density) && !isnan(row.C_overlap_density), df)

    x = valid_data.mean_gap_size
    y = valid_data.C_overlap_density

    # Objective function to minimize (SSE)
    # params = [base_rate, slope, threshold]
    function loss(params)
        base_rate, slope, threshold = params
        # Prevent negative parameters for threshold and slope
        if slope < 0 || threshold < 0 || base_rate < 0
            return Inf
        end
        preds = base_rate .+ slope .* max.(0.0, threshold .- x)
        return sum((y .- preds).^2)
    end

    # Initial guesses:
    # base_rate: minimum overlap density
    # slope: a guess
    # threshold: median of gap sizes as a starting point
    init_base = max(0.0, minimum(y))
    init_slope = 1.0
    init_thresh = median(x)

    # Optimize
    println("Optimizing piecewise threshold model...")
    res = optimize(loss, [init_base, init_slope, init_thresh], NelderMead())

    best_params = Optim.minimizer(res)
    best_base_rate, best_slope, best_threshold = best_params

    println("\n--- Mechanistic Engine Results ---")
    println("Calibrated Base Rate: ", round(best_base_rate, digits=4))
    println("Calibrated Toxicity Slope: ", round(best_slope, digits=4))
    println("Calibrated Intergenic Spacing Threshold: ", round(best_threshold, digits=4), " bp")

    # Calculate predictions using the optimized engine
    valid_data.predicted_overlap_density = best_base_rate .+ best_slope .* max.(0.0, best_threshold .- valid_data.mean_gap_size)

    # Create an engine function
    engine_func = (gap_size) -> best_base_rate + best_slope * max(0.0, best_threshold - gap_size)

    return engine_func, best_threshold, valid_data
end

# Phase 3: CairoMakie "Evolutionary Trajectory" Plot
"""
    plot_evolutionary_trajectories(df::DataFrame, engine_func, threshold::Float64, output_file::String)

Visualizes the "Genomic Squeeze" using CairoMakie.
X-Axis: OGT, Y-Axis: Genome Size.
Background contour map indicates predicted Convergent Overlap Density.
Overlays actual data points and simulated evolutionary trajectories.
"""
function plot_evolutionary_trajectories(df::DataFrame, engine_func, threshold::Float64, output_file::String)
    println("Generating Evolutionary Trajectory visualization...")

    if !("OGT" in names(df)) || !("genome_size_mb" in names(df)) || !("mean_gap_size" in names(df))
        @error "Missing required columns (OGT, genome_size_mb, mean_gap_size) for plotting."
        return
    end

    # We want to create a contour plot where Z is predicted density.
    # To do this, we need a way to estimate mean_gap_size given (OGT, genome_size_mb).
    # A simple way for visualization purposes is a linear model for gap size, or we just
    # use the actual values and interpolate.

    # Simple model to map OGT and Genome Size to Gap Size for the background contour
    mod_gap = lm(@formula(mean_gap_size ~ genome_size_mb + OGT), df)

    x_range = range(min(minimum(df.OGT), 10.0), stop=max(maximum(df.OGT), 110.0), length=100)
    y_range = range(min(minimum(df.genome_size_mb), 0.5), stop=max(maximum(df.genome_size_mb), 10.0), length=100)

    Z = zeros(100, 100)
    for (i, x) in enumerate(x_range)
        for (j, y) in enumerate(y_range)
            # Predict gap size
            pred_gap = predict(mod_gap, DataFrame(OGT=[x], genome_size_mb=[y]))[1]

            # Predict density using mechanistic engine
            pred_density = engine_func(pred_gap)
            Z[i, j] = pred_density
        end
    end

    # Setup Figure
    fig = Figure(size = (1000, 800), fontsize = 18)
    ax = Axis(fig[1, 1],
              xlabel = "Optimal Growth Temperature (°C)",
              ylabel = "Genome Size (Mb)",
              title = "The Genomic Squeeze: Evolutionary Trajectories of Genome Reduction")

    # Background Contour
    cf = contourf!(ax, x_range, y_range, Z, levels=20, colormap=:inferno)
    Colorbar(fig[1, 2], cf, label = "Predicted C_overlap Density")

    # Plot empirical data points
    scatter!(ax, df.OGT, df.genome_size_mb, color=:black, markersize=5, alpha=0.3, label="Empirical Lineages")

    # Add simulated trajectories
    # We will trace a few theoretical lineages adapting to different niches

    # 1. Mesophile streamlining (Constant OGT, shrinking size)
    traj1_ogt = fill(37.0, 50)
    traj1_size = range(6.0, stop=1.5, length=50)
    lines!(ax, traj1_ogt, traj1_size, color=:cyan, linewidth=3, label="Trajectory 1: Mesophilic Streamlining")
    scatter!(ax, traj1_ogt[1:5:end], traj1_size[1:5:end], color=:cyan, markersize=10)

    # 2. Thermophilic adaptation & extreme streamlining
    traj2_ogt = range(40.0, stop=85.0, length=50)
    traj2_size = range(4.0, stop=1.0, length=50)
    lines!(ax, traj2_ogt, traj2_size, color=:lime, linewidth=3, label="Trajectory 2: Thermophilic Reduction")
    scatter!(ax, traj2_ogt[1:5:end], traj2_size[1:5:end], color=:lime, markersize=10)

    # 3. Psychrophilic stability
    traj3_ogt = range(25.0, stop=15.0, length=50)
    traj3_size = range(5.0, stop=3.5, length=50)
    lines!(ax, traj3_ogt, traj3_size, color=:magenta, linewidth=3, label="Trajectory 3: Psychrophilic Stability")
    scatter!(ax, traj3_ogt[1:5:end], traj3_size[1:5:end], color=:magenta, markersize=10)

    # 4. Actual Empirical Trajectory (Pick a clade, e.g. a known order with variation if possible)
    # Just draw a line connecting a few real points as a pseudo-empirical trajectory
    if nrow(df) > 5
        sorted_df = sort(df, :genome_size_mb, rev=true)
        sub_df = sorted_df[1:max(2, min(5, nrow(sorted_df))), :]
        lines!(ax, sub_df.OGT, sub_df.genome_size_mb, color=:white, linewidth=2, linestyle=:dash, label="Sample Empirical Path")
    end

    axislegend(ax, position=:rt)

    println("Saving plot to $output_file")
    save(output_file, fig)
end

end # module
