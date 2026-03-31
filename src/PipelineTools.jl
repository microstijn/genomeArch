module PipelineTools

using Taxonomy
using CSV
using DataFrames
using Statistics
using GLM
using StatsModels
using Optim

export merge_and_impute_ogt
export merge_and_impute_lifestyle
export optimize_models
export mechanistic_engine


# Phase 0: OGT Database Integration & Taxonomic Imputation
"""
    merge_and_impute_ogt(genarch_file::String, tempura_file::String, taxNodesLoc::String, taxNamesLoc::String, output_file::String)

Reads genArch output and TEMPURA OGT database. Uses strict NCBI TaxID graph traversal 
via the local taxdump to calculate clade medians, completely avoiding string matching.
Imputes missing OGT values and tracks the taxonomic distance of the imputation.
"""
function merge_and_impute_ogt(genarch_file::String, tempura_file::String, taxNodesLoc::String, taxNamesLoc::String, output_file::String)
    println("Loading GenArch and TEMPURA datasets...")
    df_genarch = CSV.File(genarch_file) |> DataFrame
    df_tempura = CSV.File(tempura_file, quoted=false, silencewarnings=true) |> DataFrame

    # Standardize TEMPURA taxId column name
    if "taxonomy_id" in names(df_tempura)
        rename!(df_tempura, :taxonomy_id => :taxId)
    end

    # Helper function to safely force IDs to integers
    function safe_taxid(x)
        ismissing(x) && return missing
        s = strip(string(x))
        val = tryparse(Int, s)
        !isnothing(val) && return val
        fval = tryparse(Float64, s)
        !isnothing(fval) && return floor(Int, fval)
        return missing
    end

    df_genarch.taxId = safe_taxid.(df_genarch.taxId)
    df_tempura.taxId = safe_taxid.(df_tempura.taxId)

    println("Loading NCBI Taxonomy database (this may take a moment)...")
    db = Taxonomy.DB(taxNodesLoc, taxNamesLoc)

    println("Mapping TEMPURA data to NCBI tree...")
    tempura_exact = Dict{Int, Float64}()          # Direct Leaf ID -> OGT
    tempura_clades = Dict{Int, Vector{Float64}}() # Parent ID -> List of OGTs

    valid_tempura_ogts = Float64[]

    # Step 1: Map every TEMPURA entry to its phylogenetic lineage
    for row in eachrow(df_tempura)
        id = row.taxId
        topt = row.Topt_ave
        if !ismissing(id) && !ismissing(topt) && !isnan(topt)
            topt_val = Float64(topt)
            tempura_exact[id] = topt_val
            push!(valid_tempura_ogts, topt_val)

            # Traverse up the tree to assign this OGT to all parent clades
            try
                tax = Taxon(id, db)
                lin = Lineage(tax)
                for rank in [:species, :genus, :family, :order, :class, :phylum, :superkingdom]
                    try
                        node_str = string(lin[rank])
                        # Taxonomy.jl formats as "1234 [rank] Name". We safely extract just the integer ID.
                        m = match(r"^(\d+)", node_str)
                        if m !== nothing
                            rank_id = parse(Int, m.captures[1])
                            if !haskey(tempura_clades, rank_id)
                                tempura_clades[rank_id] = Float64[]
                            end
                            push!(tempura_clades[rank_id], topt_val)
                        end
                    catch
                        # Rank missing from this specific lineage
                    end
                end
            catch
                # ID is obsolete or missing from local NCBI dump
            end
        end
    end

    println("Calculating taxonomic medians based on exact ID clusters...")
    clade_medians = Dict{Int, Float64}()
    for (tax_id, ogt_list) in tempura_clades
        clade_medians[tax_id] = median(ogt_list)
    end
    
    global_median = length(valid_tempura_ogts) > 0 ? median(valid_tempura_ogts) : NaN

    println("Traversing genArch genomes and imputing OGT...")
    imputed_topt = Float64[]
    imputation_level = String[]
    imputation_distance = Int[]
    
    ranks_to_check = [:species, :genus, :family, :order, :class, :phylum, :superkingdom]

    for row in eachrow(df_genarch)
        id = row.taxId
        
        # 1. Missing or invalid ID
        if ismissing(id)
            push!(imputed_topt, global_median)
            push!(imputation_level, "global")
            push!(imputation_distance, 7)
            continue
        end

        # 2. Perfect match in TEMPURA database
        if haskey(tempura_exact, id)
            push!(imputed_topt, tempura_exact[id])
            push!(imputation_level, "direct")
            push!(imputation_distance, 0)
            continue
        end

        # 3. Traverse lineage upwards to find closest phylogenetic median
        val = missing
        lvl_name = "global"
        dist_score = 7

        try
            tax = Taxon(id, db)
            lin = Lineage(tax)
            
            for (i, rank) in enumerate(ranks_to_check)
                try
                    node_str = string(lin[rank])
                    m = match(r"^(\d+)", node_str)
                    if m !== nothing
                        rank_id = parse(Int, m.captures[1])
                        if haskey(clade_medians, rank_id)
                            val = clade_medians[rank_id]
                            lvl_name = string(rank)
                            dist_score = i # 1 = species, 2 = genus, etc.
                            break
                        end
                    end
                catch
                    # Rank missing from lineage
                end
            end
        catch
            # ID not found in NCBI dump
        end

        # 4. Fallback if the entire lineage failed
        if ismissing(val)
            push!(imputed_topt, global_median)
            push!(imputation_level, "global")
            push!(imputation_distance, 7)
        else
            push!(imputed_topt, val)
            push!(imputation_level, lvl_name)
            push!(imputation_distance, dist_score)
        end
    end

    df_genarch.OGT = imputed_topt
    df_genarch.imputation_level = imputation_level
    df_genarch.imputation_distance = imputation_distance

    println("Writing mapped and imputed data to $output_file...")
    CSV.write(output_file, df_genarch)
    
    return df_genarch
end

# Phase 0.5: Lifestyle Metadata Integration & Imputation
"""
    merge_and_impute_lifestyle(df::DataFrame, env_file::String, taxNodesLoc::String, taxNamesLoc::String)

Reads TSV environment data, assigns binary lifestyle scores (1.0 = Free, 0.0 = Host), 
and imputes missing values using a strict NCBI TaxID graph traversal. 
Tracks the taxonomic distance of the imputation.
"""
function merge_and_impute_lifestyle(df::DataFrame, env_file::String, taxNodesLoc::String, taxNamesLoc::String)
    println("Loading environment metadata...")
    df_env = CSV.File(env_file, delim='\t') |> DataFrame
    
    # Group by accession and collapse multiple environments into a single string
    env_grouped = combine(groupby(df_env, :accession), :environment => (x -> join(skipmissing(x), " | ")) => :environment)
    
    # Map to strict 1.0 (free) or 0.0 (host) based on regex
    host_regex = r"host|human|rumen|feces|gut|symbiont|pathogen|sheep|cow|pig|animal|plant|blood"i
    env_grouped.known_lifestyle = [
        ismissing(env) ? missing :
        occursin(host_regex, env) ? 0.0 : 1.0 
        for env in env_grouped.environment
    ]
    
    println("Merging environment data via Accession...")
    df_merged = leftjoin(df, env_grouped, on=:accession, makeunique=true)
    
    println("Loading NCBI Taxonomy database for Lifestyle Imputation...")
    db = Taxonomy.DB(taxNodesLoc, taxNamesLoc)
    
    # Map known lifestyles to clades
    clade_lifestyles = Dict{Int, Vector{Float64}}()
    valid_known = Float64[]
    
    for row in eachrow(df_merged)
        if !ismissing(row.known_lifestyle) && !ismissing(row.taxId)
            val = Float64(row.known_lifestyle)
            push!(valid_known, val)
            try
                tax = Taxon(row.taxId, db)
                lin = Lineage(tax)
                for rank in [:species, :genus, :family, :order, :class, :phylum, :superkingdom]
                    try
                        node_str = string(lin[rank])
                        m = match(r"^(\d+)", node_str)
                        if m !== nothing
                            rank_id = parse(Int, m.captures[1])
                            if !haskey(clade_lifestyles, rank_id)
                                clade_lifestyles[rank_id] = Float64[]
                            end
                            push!(clade_lifestyles[rank_id], val)
                        end
                    catch
                    end
                end
            catch
            end
        end
    end
    
    # Calculate clade means (This creates the "Probability of being free-living")
    println("Calculating taxonomic lifestyle probabilities...")
    clade_means = Dict{Int, Float64}()
    for (tax_id, vals) in clade_lifestyles
        clade_means[tax_id] = mean(vals)
    end
    
    global_mean = length(valid_known) > 0 ? mean(valid_known) : 0.5
    
    println("Traversing genArch genomes and imputing Lifestyle...")
    imputed_lifestyle = Float64[]
    imp_level = String[]
    imp_dist = Int[]
    
    ranks_to_check = [:species, :genus, :family, :order, :class, :phylum, :superkingdom]
    
    for row in eachrow(df_merged)
        # 1. We have direct empirical data
        if !ismissing(row.known_lifestyle)
            push!(imputed_lifestyle, Float64(row.known_lifestyle))
            push!(imp_level, "direct")
            push!(imp_dist, 0)
            continue
        end
        
        # 2. Missing TaxID
        if ismissing(row.taxId)
            push!(imputed_lifestyle, global_mean)
            push!(imp_level, "global")
            push!(imp_dist, 7)
            continue
        end
        
        # 3. Traverse lineage upwards to find closest phylogenetic mean
        val = missing
        lvl_name = "global"
        dist_score = 7
        
        try
            tax = Taxon(row.taxId, db)
            lin = Lineage(tax)
            
            for (i, rank) in enumerate(ranks_to_check)
                try
                    node_str = string(lin[rank])
                    m = match(r"^(\d+)", node_str)
                    if m !== nothing
                        rank_id = parse(Int, m.captures[1])
                        if haskey(clade_means, rank_id)
                            val = clade_means[rank_id]
                            lvl_name = string(rank)
                            dist_score = i
                            break
                        end
                    end
                catch
                end
            end
        catch
        end
        
        # 4. Fallback
        if ismissing(val)
            push!(imputed_lifestyle, global_mean)
            push!(imp_level, "global")
            push!(imp_dist, 7)
        else
            push!(imputed_lifestyle, val)
            push!(imp_level, lvl_name)
            push!(imp_dist, dist_score)
        end
    end
    
    df_merged.is_free_living = imputed_lifestyle
    df_merged.lifestyle_imputation_level = imp_level
    df_merged.lifestyle_imputation_distance = imp_dist
    
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

    # --- NEW: Imputation Weighting Engine ---
    println("Calculating composite imputation weights...")
    if "imputation_distance" in names(df_clean) && "lifestyle_imputation_distance" in names(df_clean)
        # Calculate inverse distance weights for both variables
        w_ogt = 1.0 ./ (1.0 .+ df_clean.imputation_distance)
        w_life = 1.0 ./ (1.0 .+ df_clean.lifestyle_imputation_distance)
        
        # Multiply them to get a composite confidence score
        raw_weights = w_ogt .* w_life
        
        # Normalize weights so they sum to the actual sample size (N).
        # This is strictly required so the log-likelihood and AICc math remain valid.
        wts = raw_weights .* (nrow(df_clean) / sum(raw_weights))
    else
        @warn "Imputation distance columns missing. Defaulting to unweighted OLS."
        wts = ones(nrow(df_clean))
    end
    df_clean.model_weights = wts

    println("Fitting Model 1 (The Rigid Baseline)...")
    # k = 2 (Intercept + genome_size_mb)
    mod1 = lm(@formula(C_overlap_density ~ genome_size_mb), df_clean, wts=wts)
    
    println("Fitting Model 2 (The Thermal Hypothesis)...")
    # k = 3 (Intercept + genome_size_mb + OGT)
    mod2 = lm(@formula(C_overlap_density ~ genome_size_mb + OGT), df_clean, wts=wts)
    
    println("Fitting Model 3 (The Decoupled Lifestyle)...")
    # k = 4 (Intercept + genome_size_mb + OGT + genome_size_mb & is_free_living)
    mod3 = lm(@formula(C_overlap_density ~ genome_size_mb + OGT + genome_size_mb & is_free_living), df_clean, wts=wts)


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
        println("Column 'mean_gap_size' not found. Calculating it dynamically...")
        if "genome_size" in names(df) && "p_gene_length_sum" in names(df) && "n_gene_length_sum" in names(df) && "total_genes" in names(df)
            # Total sequence length occupied by genes
            total_gene_len = df.p_gene_length_sum .+ df.n_gene_length_sum
            
            # The remaining non-coding space, divided by the number of gaps (which scales with total_genes)
            df.mean_gap_size = (df.genome_size .- total_gene_len) ./ df.total_genes
        else
            @error "Missing required length columns to calculate 'mean_gap_size'."
            return nothing
        end
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

end # module
