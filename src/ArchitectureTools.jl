# src/ArchitectureTools.jl

module ArchitectureTools

export calculate_architecture

# Required packages
using CSV
using DataFrames
using GFF3
using Glob # efficient data localisation
using GenomicFeatures
using ProgressMeter 
using Base.Threads # multithread support
using IntervalTrees # need for efficient overlap detection
using Statistics 

#-----------------------------------------------------------------
# Core Helper Functions
#-----------------------------------------------------------------

function length_interval(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    return array_end .- array_start
end

function space_between(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    if length(array_start) < 2
        return Int[]
    end
    return array_start[2:end] .- array_end[1:end-1]
end

function count_overlaps(length_array::AbstractArray{T}) where T<:Int
    uni_overlaps = @view length_array[length_array .<= 0]
    total_uni_overlaps = length(uni_overlaps)
    # Return the absolute sum for a positive physical length
    sum_uni_overlaps = total_uni_overlaps == 0 ? 0 : sum(abs.(uni_overlaps))
    return total_uni_overlaps, sum_uni_overlaps
end

function get_gff_file_paths(gff_dir::String)
    manifest_path = joinpath(gff_dir, "gff_manifest.txt")
    if isfile(manifest_path)
        println("Reading GFF file list from existing manifest: $manifest_path")
        return readlines(manifest_path)
    else
        println("Manifest not found. Searching for GFF files to create one (this might be slow)...")
        file_paths = glob("**/*.gff*", gff_dir)
        if isempty(file_paths)
            @warn "No GFF files found to create a manifest."
            return String[]
        end
        println("Found $(length(file_paths)) files. Writing to manifest for future runs...")
        try
            open(manifest_path, "w") do io
                for path in file_paths
                    println(io, path)
                end
            end
        catch e
            @error "Could not write manifest file at $manifest_path. Proceeding without it." exception=(e, catch_backtrace())
        end
        return file_paths
    end
end

function bidirectional_overlaps(set1, set2)
    overlap_vect_divergent = Int[]
    overlap_vect_convergent = Int[]

    for overlap in eachoverlap(set1, set2)
        interval_neg = overlap[1] # From set1 (negative strand)
        interval_pos = overlap[2] # From set2 (positive strand)

        start_neg, end_neg = leftposition(interval_neg), rightposition(interval_neg)
        start_pos, end_pos = leftposition(interval_pos), rightposition(interval_pos)

        if end_pos > start_neg && start_pos < start_neg
             push!(overlap_vect_convergent, end_pos - start_neg)
        elseif end_neg > start_pos && start_neg < start_pos
             push!(overlap_vect_divergent, end_neg - start_pos)
        end
    end

    len_con = length(overlap_vect_convergent)
    len_di = length(overlap_vect_divergent)
    sum_vect_con = isempty(overlap_vect_convergent) ? 0 : sum(overlap_vect_convergent)
    sum_vect_di = isempty(overlap_vect_divergent) ? 0 : sum(overlap_vect_divergent)

    return len_con, len_di, sum_vect_con, sum_vect_di
end


#-----------------------------------------------------------------
# Advanced Architectural Analysis Functions
#-----------------------------------------------------------------

function analyze_intergenic_distances(gaps::AbstractArray{T}) where T<:Int
    spacers = filter(x -> x > 0, gaps)
    if isempty(spacers)
        return (mean=0.0, median=0.0, std=0.0)
    end
    return (
        mean=mean(spacers),
        median=median(spacers),
        std=std(spacers)
    )
end

function analyze_operons(genes::Vector{GFF3.Record}; max_operon_gap=200)
    operons = []
    current_operon = []
    
    sorted_genes = sort(genes, by=GFF3.seqstart)

    for i in eachindex(sorted_genes)
        if isempty(current_operon)
            push!(current_operon, sorted_genes[i])
        else
            last_gene = current_operon[end]
            current_gene = sorted_genes[i]
            
            if GFF3.strand(last_gene) == GFF3.strand(current_gene) &&
               (GFF3.seqstart(current_gene) - GFF3.seqend(last_gene)) <= max_operon_gap
                push!(current_operon, current_gene)
            else
                if length(current_operon) > 1
                    push!(operons, current_operon)
                end
                current_operon = [current_gene]
            end
        end
    end
    if length(current_operon) > 1
        push!(operons, current_operon)
    end
    
    genes_in_operons = sum(length, operons; init=0)
    operonicity_score = length(genes) > 0 ? genes_in_operons / length(genes) * 100 : 0.0
    mean_operon_size = !isempty(operons) ? genes_in_operons / length(operons) : 0.0
    
    return (
        operon_nr = length(operons),
        operonicity_score = operonicity_score,
        mean_operon_size = mean_operon_size
    )
end

function analyze_local_arrangements(sorted_genes; max_spacing=500)
    divergent_pairs = 0
    convergent_pairs = 0

    for i in 1:(length(sorted_genes) - 1)
        gene1 = sorted_genes[i]
        gene2 = sorted_genes[i+1]
        
        spacing = GFF3.seqstart(gene2) - GFF3.seqend(gene1)
        if 0 < spacing <= max_spacing
            if GFF3.strand(gene1) == STRAND_NEG && GFF3.strand(gene2) == STRAND_POS
                divergent_pairs += 1
            elseif GFF3.strand(gene1) == STRAND_POS && GFF3.strand(gene2) == STRAND_NEG
                convergent_pairs += 1
            end
        end
    end
    return (divergent_pairs_nr=divergent_pairs, convergent_pairs_nr=convergent_pairs)
end

function calculate_density_gradient(genes, contig_size; window_size=100000)
    if contig_size < window_size || isempty(genes)
        return 0.0
    end
    
    densities = Float64[]
    for start in 1:window_size:contig_size
        stop = min(start + window_size - 1, contig_size)
        actual_window_size = stop - start + 1
        genes_in_window = filter(g -> GFF3.seqstart(g) >= start && GFF3.seqend(g) <= stop, genes)
        density = (length(genes_in_window) / actual_window_size) * 1000.0
        push!(densities, density)
    end
    
    return isempty(densities) ? 0.0 : std(densities)
end

#-----------------------------------------------------------------
# Main Exported Function
#-----------------------------------------------------------------
function calculate_architecture(gff_dir::String, output_file::String)
    if !isdir(gff_dir)
        @error "Input directory not found: $gff_dir"
        return
    end

    file_paths = get_gff_file_paths(gff_dir)
    n_files = length(file_paths)
    if isempty(file_paths)
        @warn "No GFF files to process. Exiting."
        return
    end
    println("Found $n_files GFF files. Starting analysis on $(nthreads()) threads...")

    thread_results = [DataFrame() for _ in 1:nthreads()]
    p = Progress(n_files, "Processing GFF files...")

    @threads for file_path in file_paths
        thread_id = threadid()
        genome_name = basename(dirname(file_path))

        try
            contig_features = Dict{String, Vector{GFF3.Record}}()
            contig_sizes = Dict{String, Int}()

            reader = GFF3.Reader(open(file_path, "r"))
            for record in reader
                seqid = GFF3.seqid(record)
                if GFF3.featuretype(record) == "region"
                    contig_sizes[seqid] = GFF3.seqend(record)
                elseif GFF3.featuretype(record) == "gene"
                    get!(() -> GFF3.Record[], contig_features, seqid)
                    push!(contig_features[seqid], record)
                end
            end
            close(reader)

            for (contig_id, contig_size) in contig_sizes
                genes = get(contig_features, contig_id, GFF3.Record[])
                all_genes_sorted = sort(genes, by=GFF3.seqstart)

                p_genes = filter(g -> GFF3.strand(g) == STRAND_POS, all_genes_sorted)
                p_starts = GFF3.seqstart.(p_genes)
                p_ends = GFF3.seqend.(p_genes)
                p_gaps = space_between(p_starts, p_ends)
                p_gene_nr = length(p_genes)
                p_gene_length_sum = isempty(p_genes) ? 0 : sum(length_interval(p_starts, p_ends))

                n_genes = filter(g -> GFF3.strand(g) == STRAND_NEG, all_genes_sorted)
                n_starts = GFF3.seqstart.(n_genes)
                n_ends = GFF3.seqend.(n_genes)
                n_gaps = space_between(n_starts, n_ends)
                n_gene_nr = length(n_genes)
                n_gene_length_sum = isempty(n_genes) ? 0 : sum(length_interval(n_starts, n_ends))

                p_U_overlap_nr, p_U_overlap_length_sum = isempty(p_gaps) ? (0, 0) : count_overlaps(p_gaps)
                n_U_overlap_nr, n_U_overlap_length_sum = isempty(n_gaps) ? (0, 0) : count_overlaps(n_gaps)
                p_gap_length_sum = isempty(p_gaps) ? 0 : sum(filter(x -> x > 0, p_gaps))
                n_gap_length_sum = isempty(n_gaps) ? 0 : sum(filter(x -> x > 0, n_gaps))
                
                if !isempty(p_genes) && !isempty(n_genes)
                    p_intervals = GenomicFeatures.Interval.(Ref(contig_id), p_starts, p_ends, Ref(STRAND_POS), 1:p_gene_nr)
                    n_intervals = GenomicFeatures.Interval.(Ref(contig_id), n_starts, n_ends, Ref(STRAND_NEG), 1:n_gene_nr)
                    # Changed sort=true to the positional argument 'true'
                    C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum = bidirectional_overlaps(
                        IntervalCollection(n_intervals, true), 
                        IntervalCollection(p_intervals, true)
                    )
                else
                    C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum = 0, 0, 0, 0
                end
                
                p_gap_stats = analyze_intergenic_distances(p_gaps)
                n_gap_stats = analyze_intergenic_distances(n_gaps)
                operon_stats = analyze_operons(all_genes_sorted)
                strand_asymmetry = isempty(genes) ? 0.5 : p_gene_nr / length(genes)
                local_arrangements = analyze_local_arrangements(all_genes_sorted)
                density_gradient_std = calculate_density_gradient(all_genes_sorted, contig_size)

                push!(thread_results[thread_id], (
                    genome_name=genome_name, contig_name=contig_id, contig_size=contig_size,
                    p_gene_nr=p_gene_nr, p_gene_length_sum=p_gene_length_sum,
                    n_gene_nr=n_gene_nr, n_gene_length_sum=n_gene_length_sum,
                    p_U_overlap_nr=p_U_overlap_nr, p_U_overlap_length_sum=p_U_overlap_length_sum,
                    n_U_overlap_nr=n_U_overlap_nr, n_U_overlap_length_sum=n_U_overlap_length_sum,
                    p_gap_length_sum=p_gap_length_sum, n_gap_length_sum=n_gap_length_sum,
                    C_overlap_nr=C_overlap_nr, D_overlap_nr=D_overlap_nr,
                    C_length_sum=C_length_sum, D_length_sum=D_length_sum,
                    p_gap_mean=p_gap_stats.mean, p_gap_median=p_gap_stats.median, p_gap_std=p_gap_stats.std,
                    n_gap_mean=n_gap_stats.mean, n_gap_median=n_gap_stats.median, n_gap_std=n_gap_stats.std,
                    operon_nr=operon_stats.operon_nr, operonicity_score=operon_stats.operonicity_score,
                    mean_operon_size=operon_stats.mean_operon_size, strand_asymmetry=strand_asymmetry,
                    divergent_pairs_nr=local_arrangements.divergent_pairs_nr,
                    convergent_pairs_nr=local_arrangements.convergent_pairs_nr,
                    gene_density_gradient_std=density_gradient_std
                ), cols=:union)
            end
        catch e
             @error "Failed to process file: $file_path" exception=(e, catch_backtrace())
        end
        next!(p)
    end

    println("\nProcessing complete. Consolidating results...")
    if all(isempty, thread_results)
        @warn "No results were generated. The output file will be empty."
        return
    end
    final_dataframe = vcat(thread_results...)

    mkpath(dirname(output_file))
    println("Writing results to: $output_file")
    try
        CSV.write(output_file, final_dataframe)
    catch e
        @error "Failed to write CSV file at $output_file." exception=(e, catch_backtrace())
    end

    println("Architecture calculation complete.")
end

end # end of module