module AnalysisTools

using DataFrames
using CSV
using Statistics

export consolidate_to_genomes

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

end # end of module

