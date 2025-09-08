module AnalysisTools

using DataFrames
using CSV
using Statistics
using Phylo
using GLM 
using StatsModels
#using AbstractTrees

export consolidate_to_genomes
export calculate_density_score
export perform_pic_analysis

"""
    build_tree_from_lineages(lineages; branchlength=1.0)

Constructs a Phylo.jl tree object from a vector of taxonomic lineages.
Each lineage should be a tuple of strings (e.g., from superkingdom to species).
"""
function build_tree_from_lineages(lineages; branchlength = 1.0)
    tree = NamedTree()
    addnode!(tree, "root")

    for lineage in lineages
        parent_node = "root"
        path = ""
        
        # Filter out both missing values and the string "NA" from the lineage
        filtered_lineage = filter(x -> !ismissing(x) && String(x) != "NA", lineage)

        # Skip this entry if the lineage is empty after filtering
        if isempty(filtered_lineage)
            continue
        end
        
        internal_lineage = filtered_lineage[1:end-1]
        tip_name = filtered_lineage[end]

        # Build the internal node path for this lineage
        for taxon_raw in internal_lineage
            taxon = String(taxon_raw)
            current_path = isempty(path) ? taxon : path * "|" * taxon
            if !hasnode(tree, current_path)
                addnode!(tree, current_path)
                addbranch!(tree, parent_node, current_path, branchlength)
            end
            parent_node = current_path
            path = current_path
        end
        
        # Add the final tip (leaf node) to the tree
        if !ismissing(tip_name)
            tip_name_str = String(tip_name)
            if !hasnode(tree, tip_name_str)
                addnode!(tree, tip_name_str)
                addbranch!(tree, parent_node, tip_name_str, branchlength)
            end
        end
    end
    return tree
end

"""
    calculate_pic(tree::AbstractTree, tip_values::Dict{String, Float64})

Calculates Felsenstein's Phylogenetic Independent Contrasts (PICs). This
implementation is robust to non-bifurcating trees (polytomies and unary nodes).
"""
function calculate_pic(tree::AbstractTree, tip_values::Dict{String, Float64})
    contrasts = Float64[]
    node_values = Dict{String, Float64}()
    
    # This is a recursive helper function that performs a post-order traversal.
    function postorder_traverse(nodename)
        # Base case: if it's a leaf, its value is the trait value from the input data.
        # CORRECTED: Explicitly qualify all tree functions with `Phylo.` to avoid conflicts.
        if Phylo.isleaf(tree, nodename)
            node_values[nodename] = get(tip_values, nodename, 0.0)
            parent = Phylo.getparent(tree, nodename)
            return Phylo.getlength(tree, Phylo.getbranch(tree, parent, nodename))
        end

        # Recursive step: process all children first.
        child_vals = Float64[]
        child_bls = Float64[]
        for child in Phylo.getchildren(tree, nodename)
            push!(child_bls, postorder_traverse(child))
            push!(child_vals, node_values[child])
        end

        # After all children are processed, calculate the value for the current internal node.
        if length(child_vals) == 2
            v1, v2 = child_vals
            bl1, bl2 = child_bls
            denominator = sqrt(bl1 + bl2)
            if denominator > 0
                push!(contrasts, (v1 - v2) / denominator)
            end
            node_values[nodename] = ((v1 * bl2) + (v2 * bl1)) / (bl1 + bl2)
            parent_bl = !Phylo.isroot(tree, nodename) ? Phylo.getlength(tree, Phylo.getbranch(tree, Phylo.getparent(tree, nodename), nodename)) : 0.0
            return parent_bl + (bl1 * bl2) / (bl1 + bl2)
        elseif length(child_vals) == 1
            node_values[nodename] = child_vals[1]
            parent_bl = !Phylo.isroot(tree, nodename) ? Phylo.getlength(tree, Phylo.getbranch(tree, Phylo.getparent(tree, nodename), nodename)) : 0.0
            return parent_bl + child_bls[1]
        elseif length(child_vals) > 2
             #@warn "Polytomy detected at node $nodename, averaging child values."
             node_values[nodename] = mean(child_vals)
             parent_bl = !Phylo.isroot(tree, nodename) ? Phylo.getlength(tree, Phylo.getbranch(tree, Phylo.getparent(tree, nodename), nodename)) : 0.0
             return parent_bl + mean(child_bls)
        end
        return 0.0
    end

    # Start the traversal from the root
    postorder_traverse(Phylo.getnodename(tree, Phylo.getroot(tree)))
    
    return contrasts
end


"""
    consolidate_to_genomes(input_file::String, output_file::String)

Aggregates a per-contig genome architecture file to a per-genome summary.
"""
function consolidate_to_genomes(input_file::String, output_file::String)
    df = CSV.File(input_file) |> DataFrame

    # Define columns for simple summation and for weighted averaging
    sum_cols = [
        :contig_size, :p_gene_nr, :n_gene_nr, :p_gene_length_sum, :n_gene_length_sum,
        :p_U_overlap_nr, :p_U_overlap_length_sum, :n_U_overlap_nr, :n_U_overlap_length_sum,
        :C_overlap_nr, :C_length_sum, :D_overlap_nr, :D_length_sum,
        :p_gap_length_sum, :n_gap_length_sum, :operon_nr,
        :divergent_pairs_nr, :convergent_pairs_nr
    ]
    
    # Columns to be averaged, weighted by contig_size
    weighted_avg_cols = [
        :p_gap_mean, :p_gap_median, :p_gap_std,
        :n_gap_mean, :n_gap_median, :n_gap_std,
        :gene_density_gradient_std
    ]

    gdf = groupby(df, :genome_name)

    # Perform aggregation with explicit column naming to avoid conflicts
    consolidated_df = combine(gdf) do sub_df
        # --- Create a NamedTuple for all the simple sums ---
        # CORRECTED: Added parentheses around the entire generator expression
        sums = (; ((Symbol(col) => sum(sub_df[!, col])) for col in sum_cols)...)
        
        # --- Calculate weighted averages ---
        total_size = sums.contig_size
        weighted_avgs = NamedTuple()
        if total_size > 0
            # CORRECTED: Added parentheses around the entire generator expression
            weighted_avgs = (; ((Symbol(col) => sum(sub_df[!, col] .* sub_df.contig_size) / total_size) for col in weighted_avg_cols)...)
        else
            # Default to 0 if total size is 0 to avoid division by zero
            # CORRECTED: Added parentheses around the entire generator expression
            weighted_avgs = (; ((Symbol(col) => 0.0) for col in weighted_avg_cols)...)
        end
        
        # --- Recalculate ratios at the genome level ---
        total_genes = sums.p_gene_nr + sums.n_gene_nr
        
        recalculated_ratios = (
            strand_asymmetry = total_genes > 0 ? sums.p_gene_nr / total_genes : 0.5,
            # Operonicity is (total genes in operons) / (total genes)
            operonicity_score = total_genes > 0 ? (sum(sub_df.mean_operon_size .* sub_df.operon_nr) / total_genes) * 100 : 0.0,
            # Mean operon size is a weighted average of the mean sizes from each contig
            mean_operon_size = sums.operon_nr > 0 ? sum(sub_df.mean_operon_size .* sub_df.operon_nr) / sums.operon_nr : 0.0
        )
        
        # Merge all results into a single NamedTuple for the new row
        other_sums = Base.structdiff(sums, (contig_size = nothing,))
        return merge((genome_size=sums.contig_size,), other_sums, weighted_avgs, recalculated_ratios)
    end
    
    CSV.write(output_file, consolidated_df)
    println("Consolidation complete. Per-genome results written to: $output_file")
end

# --- New Function ---
"""
    calculate_density_score(input_file::String, output_file::String)

Reads a consolidated genome metrics file, engineers features, and calculates
a composite Density Score based on normalized architectural metrics.
"""
function calculate_density_score(input_file::String, output_file::String)
    df = CSV.File(input_file) |> DataFrame
    println("Engineering features for density score calculation...")
    
    # --- 1. Feature Engineering (as described in the paper) ---
    df.total_genes = df.p_gene_nr .+ df.n_gene_nr
    df = filter(row -> row.total_genes > 0, df) # Avoid division by zero

    overlap_cols = [:p_U_overlap_nr, :n_U_overlap_nr, :C_overlap_nr, :D_overlap_nr]
    df.total_overlaps = sum.(eachrow(df[!, overlap_cols]))
    df.overlaps_per_1000_genes = (df.total_overlaps ./ df.total_genes) .* 1000

    df.mean_gap_size = zeros(nrow(df))
    for r in eachrow(df)
        gap_sum = r.p_gap_length_sum + r.n_gap_length_sum
        gene_count_for_gaps = r.total_genes
        if gene_count_for_gaps > 2
            r.mean_gap_size = gap_sum / (gene_count_for_gaps - 2)
        end
    end

    # --- 2. Standardization (Z-score) ---
    metrics_for_score = [:overlaps_per_1000_genes, :mean_gap_size, :operonicity_score]
    
    for metric in metrics_for_score
        # Impute any potential non-finite values with the column median
        if any(!isfinite, df[!, metric])
            m = median(filter(isfinite, df[!, metric]))
            df[!, metric] = [isfinite(x) ? x : m for x in df[!, metric]]
        end
        
        μ = mean(df[!, metric])
        σ = std(df[!, metric])
        z_score_col = Symbol("z_", metric)
        df[!, z_score_col] = σ > 0 ? (df[!, metric] .- μ) ./ σ : 0.0
    end

    # --- 3. Calculate Composite Density Score ---
    println("Calculating composite Density Score...")
    df.DensityScore = df.z_overlaps_per_1000_genes .- df.z_mean_gap_size .+ df.z_operonicity_score

    # --- 4. Save Results ---
    # We keep the engineered features but remove the intermediate z-score columns
    final_df = select(df, Not(r"^z_"))
    CSV.write(output_file, final_df)
    println("Density score calculation complete. Results written to: $output_file")
    
    return final_df
end

"""
    perform_pic_analysis(data_file::String, taxonomy_file::String, environment_file::String)

Performs a Phylogenetic Independent Contrasts (PIC) analysis to test the
correlation between DensityScore and Environment while correcting for phylogeny.
"""
function perform_pic_analysis(data_file::String, taxonomy_file::String, environment_file::String)
    println("Starting Phylogenetic Independent Contrasts (PIC) analysis...")
    
    # --- 1. Load All Data Sources ---
    df = CSV.File(data_file) |> DataFrame
    tax_df = CSV.File(taxonomy_file) |> DataFrame
    env_df = CSV.File(environment_file, delim = '\t') |> DataFrame

    # --- 2. Standardize Column Names for Joins ---
    tax_df_names = names(tax_df)
    if "genome_name" in tax_df_names
    elseif "accession" in tax_df_names
        rename!(tax_df, :accession => :genome_name)
    elseif "assembly_accession" in tax_df_names
        rename!(tax_df, :assembly_accession => :genome_name)
    else
        error("Taxonomy file must contain 'genome_name', 'accession', or 'assembly_accession'.")
    end

    if "accession" in names(env_df)
        rename!(env_df, :accession => :genome_name)
    end
    
    # --- 3. Process and Classify Environments ---
    host_associated_genomes = Set(filter(row -> !ismissing(row.environment) && occursin(r"host|human|rumen|feces|gut"i, row.environment), env_df).genome_name)
    env_classification = DataFrame(genome_name = unique(tax_df.genome_name))
    env_classification.is_host_associated = [name in host_associated_genomes ? 1.0 : 0.0 for name in env_classification.genome_name]

    # --- 4. Merge the three DataFrames ---
    merged_df = innerjoin(df, tax_df, on = :genome_name)
    merged_df = leftjoin(merged_df, env_classification, on = :genome_name)
    merged_df.is_host_associated = coalesce.(merged_df.is_host_associated, 0.0)

    # --- 5. Build a Phylogenetic Tree from Taxonomy ---
    println("Constructing taxonomic tree...")
    lineages = [(r.superkingdom, r.phylum, r.class, r.order, r.family, r.genus, r.species, r.genome_name) for r in eachrow(merged_df)]
    tree = build_tree_from_lineages(lineages)
    
    tip_names = getleafnames(tree)
    
    # --- 6. Map Traits to the Tree ---
    tip_data = DataFrame(genome_name = tip_names)
    tip_data = leftjoin(tip_data, merged_df, on = :genome_name)
    
    density_scores = coalesce.(tip_data.DensityScore, median(skipmissing(tip_data.DensityScore)))
    environments = coalesce.(tip_data.is_host_associated, 0.0)

    # Create dictionaries mapping tip names to values
    density_dict = Dict(zip(tip_data.genome_name, density_scores))
    env_dict = Dict(zip(tip_data.genome_name, environments))

    # --- 7. Calculate Independent Contrasts ---
    println("Calculating independent contrasts...")
    contrasts_density = calculate_pic(tree, density_dict)
    contrasts_environment = calculate_pic(tree, env_dict)
    
    # --- 8. Perform Statistical Test ---
    println("Performing regression on contrasts...")
    # Ensure there's a 1-to-1 match between contrasts before creating the DataFrame
    min_len = min(length(contrasts_density), length(contrasts_environment))
    contrast_df = DataFrame(Density = contrasts_density[1:min_len], Environment = contrasts_environment[1:min_len])
    
    model = lm(@formula(Density ~ 0 + Environment), contrast_df)
    
    # --- 9. Report Results ---
    println("\n--- PIC Analysis Results ---")
    println(model)
    r_squared = r2(model)
    p_value = coeftable(model).cols[4][1]
    
    println("\nR-squared of contrasts: ", round(r_squared, digits = 4))
    println("P-value of correlation: ", p_value)
    
    if p_value < 0.05
        println("\nConclusion: The correlation between Density Score and Environment is statistically significant even after correcting for phylogenetic non-independence.")
    else
        println("\nConclusion: The correlation between Density Score and Environment is NOT statistically significant after correcting for phylogenetic non-independence.")
    end
    println("--------------------------")

    return model
end


end # end of module

