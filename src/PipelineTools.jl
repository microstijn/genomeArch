module PipelineTools

using Taxonomy
using CSV
using DataFrames
using Statistics
using GLM
using StatsModels
using Optim
using CategoricalArrays

export merge_and_impute_ogt
export merge_and_impute_lifestyle
export optimize_models
export mechanistic_engine
export mechanistic_engine_all


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

# Phase 0.5: Lifestyle Metadata Integration & Imputation (3-Tier Showdown)

using CSV
using DataFrames
using Statistics

"""
    merge_and_impute_lifestyle(df::DataFrame, env_file::String, taxNodesLoc::String, taxNamesLoc::String)

Reads TSV environment data, assigns 3-tier categorical lifestyle tiers 
("1_FreeLiving", "2_ExtracellularHost", "3_Endosymbiont"), and imputes missing 
values using strict NCBI TaxID graph traversal (majority voting per clade). 
Tracks the taxonomic distance of the imputation.
"""
function merge_and_impute_lifestyle(df::DataFrame, env_file::String, taxNodesLoc::String, taxNamesLoc::String)
    println("Loading environment metadata...")
    df_env = CSV.File(env_file, delim='\t') |> DataFrame
    
    # Group by accession and collapse multiple environments into a single string
    env_grouped = combine(groupby(df_env, :accession), :environment => (x -> join(skipmissing(x), " | ")) => :environment)
    
    # --- 1. THE 3-TIER REGEX MAPPING ---
    # Endosymbionts take priority. If it says "obligate endosymbiont of insect gut", it catches "endo" first.
    endo_regex = r"endosymbiont|intracellular|obligate"i
    host_regex = r"host|human|rumen|feces|gut|symbiont|pathogen|sheep|cow|pig|animal|plant|blood"i
    
    env_grouped.known_lifestyle_tier = [
        ismissing(env) ? missing :
        occursin(endo_regex, env) ? "3_Endosymbiont" :
        occursin(host_regex, env) ? "2_ExtracellularHost" : "1_FreeLiving"
        for env in env_grouped.environment
    ]
    
    println("Merging environment data via Accession...")
    df_merged = leftjoin(df, env_grouped, on=:accession, makeunique=true)
    
    println("Loading NCBI Taxonomy database for Lifestyle Imputation...")
    db = Taxonomy.DB(taxNodesLoc, taxNamesLoc)
    
    # --- 2. MAP KNOWN LIFESTYLES TO CLADES ---
    clade_lifestyles = Dict{Int, Vector{String}}()
    valid_known = String[]
    
    for row in eachrow(df_merged)
        if !ismissing(row.known_lifestyle_tier) && !ismissing(row.taxId)
            val = String(row.known_lifestyle_tier)
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
                                clade_lifestyles[rank_id] = String[]
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
    
    # --- 3. CALCULATE CLADE MAJORITY (MODE) ---
    println("Calculating taxonomic lifestyle majorities (Voting)...")
    clade_majority = Dict{Int, String}()
    
    # Helper to find the most frequent string in an array
    function get_mode(arr)
        counts = Dict{String, Int}()
        for v in arr
            counts[v] = get(counts, v, 0) + 1
        end
        best_k = ""
        best_v = -1
        for (k, v) in counts
            if v > best_v
                best_v = v
                best_k = k
            end
        end
        return best_k
    end

    for (tax_id, vals) in clade_lifestyles
        clade_majority[tax_id] = get_mode(vals)
    end
    
    global_majority = length(valid_known) > 0 ? get_mode(valid_known) : "1_FreeLiving"
    
    # --- 4. TRAVERSE AND IMPUTE ---
    println("Traversing genArch genomes and imputing Lifestyle...")
    imputed_tier = String[]
    imp_level = String[]
    imp_dist = Int[]
    
    ranks_to_check = [:species, :genus, :family, :order, :class, :phylum, :superkingdom]
    
    for row in eachrow(df_merged)
        # 1. We have direct empirical data
        if !ismissing(row.known_lifestyle_tier)
            push!(imputed_tier, String(row.known_lifestyle_tier))
            push!(imp_level, "direct")
            push!(imp_dist, 0)
            continue
        end
        
        # 2. Missing TaxID
        if ismissing(row.taxId)
            push!(imputed_tier, global_majority)
            push!(imp_level, "global")
            push!(imp_dist, 7)
            continue
        end
        
        # 3. Traverse lineage upwards to find closest phylogenetic majority
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
                        if haskey(clade_majority, rank_id)
                            val = clade_majority[rank_id]
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
            push!(imputed_tier, global_majority)
            push!(imp_level, "global")
            push!(imp_dist, 7)
        else
            push!(imputed_tier, val)
            push!(imp_level, lvl_name)
            push!(imp_dist, dist_score)
        end
    end
    
    # Assign the new categorical column
    df_merged.lifestyle_tier = imputed_tier
    df_merged.lifestyle_imputation_level = imp_level
    df_merged.lifestyle_imputation_distance = imp_dist
    
    # Create the backward-compatible binary column just in case Phase 1 needs it
    # Maps FreeLiving -> 1.0, everything else (ExtracellularHost & Endosymbiont) -> 0.0
    df_merged.is_free_living = [tier == "1_FreeLiving" ? 1.0 : 0.0 for tier in imputed_tier]
    
    return df_merged
end

# Phase 1: Statistical Modeling & Optimization (3-Tier Showdown with BIC)

using GLM
using DataFrames
using Statistics
using CategoricalArrays

"""
    optimize_models(df::DataFrame, threshold::Float64, overlap_col::Symbol)

Dynamically calculates derived metrics, engineers a quadratic compression penalty, 
applies composite inverse-distance weights, and evaluates 6 evolutionary models 
using Weighted Least Squares (WLS).

Evaluated using the rigorous Bayesian Information Criterion (BIC) to heavily penalize 
overfitting on large datasets, preventing localized noise (like Actinomycetota 
antisense networks) from distorting global evolutionary rules.
"""
function optimize_models(df::DataFrame, threshold::Float64, overlap_col::Symbol)
    println("\n=== Preparing Data for Optimization: $(overlap_col) ===")
    df_clean = copy(df) 
    
    # --- 1. DYNAMICALLY CALCULATE MISSING METRICS ---
    if !("genome_size_mb" in names(df_clean))
        df_clean.genome_size_mb = df_clean.genome_size ./ 1_000_000
    end
    
    if !("total_genes" in names(df_clean)) && "p_gene_nr" in names(df_clean)
        df_clean.total_genes = df_clean.p_gene_nr .+ df_clean.n_gene_nr
    end

    if !("mean_gap_size" in names(df_clean))
        println("Calculating 'mean_gap_size' dynamically...")
        total_gene_len = df_clean.p_gene_length_sum .+ df_clean.n_gene_length_sum
        total_overlap_len = df_clean.p_U_overlap_length_sum .+ 
                            df_clean.n_U_overlap_length_sum .+ 
                            df_clean.C_length_sum .+ 
                            df_clean.D_length_sum
        true_gap_space = df_clean.genome_size .- total_gene_len .+ total_overlap_len
        df_clean.mean_gap_size = true_gap_space ./ df_clean.total_genes
    end

    println("Calculating density for $(overlap_col) (per 1000 genes)...")
    df_clean.target_density = (df_clean[!, overlap_col] ./ df_clean.total_genes) .* 1000.0

    # --- 2. ENGINEER THE DYNAMIC HINGE VARIABLES ---
    println("Engineering Quadratic Compression Penalty (Threshold: $(round(threshold, digits=3)) bp)...")
    df_clean.compression_penalty = max.(0.0, threshold .- df_clean.mean_gap_size)
    df_clean.compression_penalty_sq = df_clean.compression_penalty .^ 2

    # --- 3. PREPARE 3-TIER CATEGORICAL LIFESTYLE ---
    println("Setting ExtracellularHost as the statistical baseline...")
    df_clean.lifestyle_tier = categorical(df_clean.lifestyle_tier)
    levels!(df_clean.lifestyle_tier, ["2_ExtracellularHost", "1_FreeLiving", "3_Endosymbiont"])

    # --- 4. FILTER OUT MISSING DATA ---
    cols_to_check = [:target_density, :genome_size_mb, :OGT, :lifestyle_tier, :compression_penalty, :total_genes, :mean_gap_size]
    dropmissing!(df_clean, cols_to_check)

    # --- 5. CALCULATE IMPUTATION WEIGHTS ---
    if "imputation_distance" in names(df_clean) && "lifestyle_imputation_distance" in names(df_clean)
        w_ogt = 1.0 ./ (1.0 .+ df_clean.imputation_distance)
        w_life = 1.0 ./ (1.0 .+ df_clean.lifestyle_imputation_distance)
        raw_weights = w_ogt .* w_life
        wts = raw_weights .* (nrow(df_clean) / sum(raw_weights)) 
    else
        wts = ones(nrow(df_clean))
    end
    df_clean.model_weights = wts

    # --- 6. HELPER FUNCTION FOR BIC (The Ruthless Overfit Penalty) ---
    function calc_bic(m)
        k = length(coef(m)) + 1 # Number of parameters + 1 for variance
        n = nobs(m)             # Sample size
        ll = loglikelihood(m)
        return k * log(n) - 2 * ll
    end

    # --- 7. FIT THE 6 MODELS (Safely Named to Avoid 'mod1' function clash) ---
    println("Fitting Model 1 (The Rigid Baseline)...")
    model1 = lm(@formula(target_density ~ genome_size_mb), df_clean, wts=wts)
    
    println("Fitting Model 2 (The Thermal Hypothesis)...")
    model2 = lm(@formula(target_density ~ genome_size_mb + OGT), df_clean, wts=wts)
    
    println("Fitting Model 3 (The Decoupled Lifestyle)...")
    model3 = lm(@formula(target_density ~ genome_size_mb + OGT + genome_size_mb & lifestyle_tier), df_clean, wts=wts)

    println("Fitting Model 4 (The Physical Threshold)...")
    model4 = lm(@formula(target_density ~ genome_size_mb + OGT + genome_size_mb & lifestyle_tier + 
                         compression_penalty), df_clean, wts=wts)

    println("Fitting Model 5 (The Ecological Squeeze)...")
    model5 = lm(@formula(target_density ~ genome_size_mb + OGT + genome_size_mb & lifestyle_tier + 
                         compression_penalty + compression_penalty & lifestyle_tier), df_clean, wts=wts)

    println("Fitting Model 6 (The 3-Tier Quadratic Squeeze)...")
    model6 = lm(@formula(target_density ~ genome_size_mb + OGT + genome_size_mb & lifestyle_tier + 
                         compression_penalty + compression_penalty_sq + 
                         compression_penalty_sq & lifestyle_tier), df_clean, wts=wts)

    # --- 8. COMPARE SCORES USING BIC ---
    bics = [calc_bic(model1), calc_bic(model2), calc_bic(model3), calc_bic(model4), calc_bic(model5), calc_bic(model6)]
    models_array = [model1, model2, model3, model4, model5, model6]
    
    println("\n--- BIC Model Comparison (Strict Complexity Penalty) ---")
    println("Model 1 (Baseline):             BIC = $(round(bics[1], digits=2))")
    println("Model 2 (Thermal Hypothesis):   BIC = $(round(bics[2], digits=2))")
    println("Model 3 (Decoupled Lifestyle):  BIC = $(round(bics[3], digits=2))")
    println("Model 4 (Physical Threshold):   BIC = $(round(bics[4], digits=2))")
    println("Model 5 (Ecological Squeeze):   BIC = $(round(bics[5], digits=2))")
    println("Model 6 (3-Tier Quad Squeeze):  BIC = $(round(bics[6], digits=2))")
    
    best_idx = argmin(bics)
    best_model = models_array[best_idx]
    
    println("\n-> The winning model under rigorous BIC is: Model $best_idx")
    println(coeftable(best_model))
    println("-> Winning Model R² = ", round(r2(best_model), digits=4))

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
        println("Calculating 'mean_gap_size' from strand-specific positive gaps...")
        
        # 1. Total number of genes per strand
        p_genes = df.p_gene_nr
        n_genes = df.n_gene_nr
        total_genes = p_genes .+ n_genes
        df.total_genes = total_genes
        
        # 2. Safely calculate the weighted global mean gap size (excluding overlaps)
        # Using the p_gap_mean and n_gap_mean which (based on your methods) 
        # should already have negative distances filtered out.
        df.mean_gap_size = ((df.p_gap_mean .* p_genes) .+ (df.n_gap_mean .* n_genes)) ./ total_genes 
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
        preds = base_rate .+ slope .* max.(0.0, threshold .- x)
        return sum((y .- preds).^2)
    end
    
    # Initial guesses
    init_base = max(0.0, minimum(y))
    init_slope = 1.0
    # Safe init: ensure the median starting point is actually inside our biological bounds
    init_thresh = clamp(median(x), 50.0, 250.0) 
    init_params = [init_base, init_slope, init_thresh]
    
    # ---------------------------------------------------------
    # THE FIX: STRICT BIOLOGICAL GUARDRAILS
    # ---------------------------------------------------------
    # [base_rate, slope, threshold]
    lower_bounds = [0.0, 0.0, 50.0]       # Threshold CANNOT go below 50 bp
    upper_bounds = [300.0, 200.0, 250.0]  # Threshold CANNOT go above 250 bp
    
    println("Optimizing piecewise threshold model with strict biological bounds...")
    res = optimize(loss, lower_bounds, upper_bounds, init_params, Fminbox(NelderMead()))
    
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



function mechanistic_engine_all(df::DataFrame)
    println("Initializing the Multi-Target Decoupled Mechanistic Engine...")
    
    # 1. Calculate Mean Gap Size if missing
    if !("mean_gap_size" in names(df))
        println("Calculating 'mean_gap_size' dynamically...")
        df.total_genes = df.p_gene_nr .+ df.n_gene_nr
        total_gene_len = df.p_gene_length_sum .+ df.n_gene_length_sum
        total_overlap_len = df.p_U_overlap_length_sum .+ df.n_U_overlap_length_sum .+ df.C_length_sum .+ df.D_length_sum
        true_gap_space = df.genome_size .- total_gene_len .+ total_overlap_len
        df.mean_gap_size = true_gap_space ./ df.total_genes
    end
    
    # 2. Calculate Densities for ALL types
    df.C_overlap_density = (df.C_overlap_nr ./ df.total_genes) .* 1000.0
    df.U_overlap_density = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
    df.D_overlap_density = (df.D_overlap_nr ./ df.total_genes) .* 1000.0
    df.Total_overlap_density = df.C_overlap_density .+ df.U_overlap_density .+ df.D_overlap_density

    # Helper function to run the optimizer on a specific target
    function optimize_target(target_col::Symbol)
        valid_data = filter(row -> !ismissing(row.mean_gap_size) && !isnan(row.mean_gap_size) && 
                                   !ismissing(row[target_col]) && !isnan(row[target_col]), df)
        
        x = Float64.(valid_data.mean_gap_size)
        y = Float64.(valid_data[!, target_col])
        
        function loss(params)
            base_rate, slope, threshold = params
            if slope < 0 || threshold < 0 || base_rate < 0
                return Inf
            end
            preds = base_rate .+ slope .* max.(0.0, threshold .- x)
            return sum((y .- preds).^2)
        end
        
        init_base = max(0.0, minimum(y))
        init_slope = 1.0
        init_thresh = median(x)
        
        res = optimize(loss, [init_base, init_slope, init_thresh], NelderMead())
        best_base_rate, best_slope, best_threshold = Optim.minimizer(res)
        
        return best_base_rate, best_slope, best_threshold
    end

    # 3. Run Optimization for each type
    targets = [:Total_overlap_density, :C_overlap_density, :U_overlap_density, :D_overlap_density]
    results = DataFrame(Overlap_Type=String[], Base_Rate=Float64[], Surge_Slope=Float64[], Threshold_bp=Float64[])
    
    for target in targets
        base, slope, thresh = optimize_target(target)
        push!(results, (string(target), base, slope, thresh))
    end
    
    println("\n=== Mechanistic Engine Results ===")
    println(results)
    
    return results
end

# Run it on your dataframe:
# results_df = mechanistic_engine_all(df_clean)

end # module
