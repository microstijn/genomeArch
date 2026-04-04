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
    return array_end .- array_start .+ 1
end

function space_between(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    if length(array_start) < 2
        return Int[]
    end
    return array_start[2:end] .- array_end[1:end-1] .- 1
end

function count_overlaps(length_array::AbstractArray{T}) where T<:Int
    uni_overlaps = @view length_array[length_array .< 0]
    abs_overlaps = abs.(uni_overlaps)
    
    total_nr = length(abs_overlaps)
    sum_len = total_nr == 0 ? 0 : sum(abs_overlaps)
    max_len = total_nr == 0 ? 0 : maximum(abs_overlaps)
    
    coupled_nr = count(x -> x == 1 || x == 4, abs_overlaps)
    deep_nr = count(x -> x > 50, abs_overlaps)
    
    in_frame_nr = count(x -> x % 3 == 0, abs_overlaps)
    out_of_frame_nr = total_nr - in_frame_nr
    abutting_nr = count(x -> x == 0, length_array)
    
    return total_nr, sum_len, max_len, coupled_nr, deep_nr, in_frame_nr, out_of_frame_nr, abutting_nr
end

function calculate_absolute_gaps(genes::Vector{GFF3.Record})
    # NOTE: Assumes `genes` is already sorted by seqstart
    if isempty(genes)
        return Int[]
    end
    
    blocks = []
    current_start = GFF3.seqstart(genes[1])
    current_end = GFF3.seqend(genes[1])
    
    for i in 2:length(genes)
        g_start = GFF3.seqstart(genes[i])
        g_end = GFF3.seqend(genes[i])
        
        if g_start <= current_end
            current_end = max(current_end, g_end)
        else
            push!(blocks, (current_start, current_end))
            current_start = g_start
            current_end = g_end
        end
    end
    push!(blocks, (current_start, current_end))
    
    gaps = Int[]
    if length(blocks) > 1
        for i in 1:(length(blocks)-1)
            gap_size = blocks[i+1][1] - blocks[i][2] - 1
            if gap_size > 0
                push!(gaps, gap_size)
            end
        end
    end
    
    return gaps
end

function analyze_overlap_chains(genes::Vector{GFF3.Record})
    # NOTE: Assumes `genes` is already sorted by seqstart
    if length(genes) < 2
        return (mean_chain_size=0.0, max_chain_size=0, total_chains=0)
    end
    
    chains = Int[]
    current_chain_length = 1
    current_end = GFF3.seqend(genes[1])
    
    for i in 2:length(genes)
        g_start = GFF3.seqstart(genes[i])
        g_end = GFF3.seqend(genes[i])
        
        if g_start <= current_end
            current_chain_length += 1
            current_end = max(current_end, g_end)
        else
            if current_chain_length > 1
                push!(chains, current_chain_length)
            end
            current_chain_length = 1
            current_end = g_end
        end
    end
    if current_chain_length > 1
        push!(chains, current_chain_length)
    end
    
    return (
        mean_chain_size = isempty(chains) ? 0.0 : mean(chains),
        max_chain_size = isempty(chains) ? 0 : maximum(chains),
        total_chains = length(chains)
    )
end

function max_overlaps_per_gene(sorted_genes::Vector{GFF3.Record})
    if isempty(sorted_genes) return 0 end
    
    # Pre-extract integers to avoid massive O(N^2) GFF parsing overhead
    starts = GFF3.seqstart.(sorted_genes)
    ends = GFF3.seqend.(sorted_genes)
    
    n = length(starts)
    max_ov = 0
    
    for i in 1:n
        ov = 0
        g1_s = starts[i]
        g1_e = ends[i]
        
        # Look ahead
        for j in (i+1):n
            if starts[j] <= g1_e
                ov += 1
            else
                break # Broken early because the array is sorted
            end
        end
        
        # Look behind
        for j in (i-1):-1:1
            if ends[j] >= g1_s
                ov += 1
            end
        end
        
        if ov > max_ov
            max_ov = ov
        end
    end
    
    return max_ov
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
        interval_neg = overlap[1]
        interval_pos = overlap[2]

        start_neg, end_neg = leftposition(interval_neg), rightposition(interval_neg)
        start_pos, end_pos = leftposition(interval_pos), rightposition(interval_pos)

        o_start = max(start_neg, start_pos)
        o_end = min(end_neg, end_pos)
        o_length = o_end - o_start + 1 
        
        if end_pos >= start_neg && start_pos < start_neg
            push!(overlap_vect_convergent, o_length)
        elseif start_pos <= end_neg && end_pos > end_neg
            push!(overlap_vect_divergent, o_length)
        elseif (start_pos >= start_neg && end_pos <= end_neg) || 
               (start_neg >= start_pos && end_neg <= end_pos)
            push!(overlap_vect_convergent, o_length)
        end
    end

    len_con = length(overlap_vect_convergent)
    len_di = length(overlap_vect_divergent)
    sum_vect_con = isempty(overlap_vect_convergent) ? 0 : sum(overlap_vect_convergent)
    sum_vect_di = isempty(overlap_vect_divergent) ? 0 : sum(overlap_vect_divergent)

    max_con = isempty(overlap_vect_convergent) ? 0 : maximum(overlap_vect_convergent)
    max_di = isempty(overlap_vect_divergent) ? 0 : maximum(overlap_vect_divergent)
    
    deep_con = count(x -> x > 50, overlap_vect_convergent)
    deep_di = count(x -> x > 50, overlap_vect_divergent)

    return len_con, len_di, sum_vect_con, sum_vect_di, max_con, max_di, deep_con, deep_di
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

function analyze_accordion_gaps(genes::Vector{GFF3.Record}; max_operon_gap=200)
    intra_gaps = Int[]
    inter_gaps = Int[]

    for strand in [STRAND_POS, STRAND_NEG]
        strand_genes = filter(g -> GFF3.strand(g) == strand, genes)
        if length(strand_genes) < 2 continue end
        
        for i in 1:(length(strand_genes)-1)
            gap = GFF3.seqstart(strand_genes[i+1]) - GFF3.seqend(strand_genes[i]) - 1
            if gap > 0
                if gap <= max_operon_gap
                    push!(intra_gaps, gap)
                else
                    push!(inter_gaps, gap)
                end
            end
        end
    end
    
    return (
        intra_gap_mean = isempty(intra_gaps) ? 0.0 : mean(intra_gaps),
        inter_gap_mean = isempty(inter_gaps) ? 0.0 : mean(inter_gaps)
    )
end

function analyze_operons(genes::Vector{GFF3.Record}; max_operon_gap=200)
    # NOTE: Assumes `genes` is already sorted by seqstart
    operons = []
    current_operon = []

    for i in eachindex(genes)
        if isempty(current_operon)
            push!(current_operon, genes[i])
        else
            last_gene = current_operon[end]
            current_gene = genes[i]
            
            if GFF3.strand(last_gene) == GFF3.strand(current_gene) &&
               (GFF3.seqstart(current_gene) - GFF3.seqend(last_gene) - 1) <= max_operon_gap
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

function analyze_local_arrangements(genes; max_spacing=500)
    # NOTE: Assumes `genes` is already sorted by seqstart
    divergent_pairs = 0
    convergent_pairs = 0

    for i in 1:(length(genes) - 1)
        gene1 = genes[i]
        gene2 = genes[i+1]
        
        spacing = GFF3.seqstart(gene2) - GFF3.seqend(gene1) - 1
        if spacing <= max_spacing
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

function calculate_coding_density(genes::Vector{GFF3.Record}, contig_size::Int)
    # NOTE: Assumes `genes` is already sorted by seqstart
    if isempty(genes) || contig_size <= 0
        return 0.0
    end
    
    current_start = GFF3.seqstart(genes[1])
    current_end = GFF3.seqend(genes[1])
    total_coding_bp = 0
    
    for i in 2:length(genes)
        g_start = GFF3.seqstart(genes[i])
        g_end = GFF3.seqend(genes[i])
        
        if g_start <= current_end
            current_end = max(current_end, g_end)
        else
            total_coding_bp += (current_end - current_start + 1)
            current_start = g_start
            current_end = g_end
        end
    end
    total_coding_bp += (current_end - current_start + 1)
    
    return (total_coding_bp / contig_size) * 100.0
end

function count_nested_genes(genes::Vector{GFF3.Record})
    if length(genes) < 2
        return 0
    end
    
    # CRITICAL: Custom sort required here. Cannot use strictly seqstart-sorted array.
    sorted_genes = sort(genes, lt=(a,b) -> begin
        sa = GFF3.seqstart(a)
        sb = GFF3.seqstart(b)
        if sa == sb
            return GFF3.seqend(a) > GFF3.seqend(b)
        else
            return sa < sb
        end
    end)
    
    nested_count = 0
    max_end_so_far = GFF3.seqend(sorted_genes[1])
    
    for i in 2:length(sorted_genes)
        g_end = GFF3.seqend(sorted_genes[i])
        
        if g_end <= max_end_so_far
            nested_count += 1
        else
            max_end_so_far = g_end
        end
    end
    
    return nested_count
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

            # Safe file reading to prevent locks
            open(file_path, "r") do io
                reader = GFF3.Reader(io)
                for record in reader
                    seqid = GFF3.seqid(record)
                    if GFF3.featuretype(record) == "region"
                        contig_sizes[seqid] = GFF3.seqend(record)
                    elseif GFF3.featuretype(record) == "gene"
                        get!(() -> GFF3.Record[], contig_features, seqid)
                        push!(contig_features[seqid], record)
                    end
                end
            end

            for (contig_id, contig_size) in contig_sizes
                genes = get(contig_features, contig_id, GFF3.Record[])
                all_genes_sorted = sort(genes, by=GFF3.seqstart)
                num_genes = length(all_genes_sorted)

                if num_genes > 0
                    all_gene_lengths = length_interval(GFF3.seqstart.(all_genes_sorted), GFF3.seqend.(all_genes_sorted))
                    mean_gene_length = mean(all_gene_lengths)
                    std_gene_length = std(all_gene_lengths)
                else
                    mean_gene_length = 0.0
                    std_gene_length = 0.0
                end
                
                strand_switches = num_genes > 1 ? count(i -> GFF3.strand(all_genes_sorted[i]) != GFF3.strand(all_genes_sorted[i+1]), 1:(num_genes-1)) : 0
                strand_switch_rate = num_genes > 0 ? strand_switches / num_genes : 0.0

                absolute_gaps = calculate_absolute_gaps(all_genes_sorted)
                absolute_gap_stats = analyze_intergenic_distances(absolute_gaps)
                
                overlap_chains = analyze_overlap_chains(all_genes_sorted)
                accordion_gaps = analyze_accordion_gaps(all_genes_sorted)
                entanglement_hubs_max = max_overlaps_per_gene(all_genes_sorted)

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

                p_U_overlap_nr, p_U_overlap_length_sum, p_U_max_len, p_U_coupled_nr, p_U_deep_nr, p_U_in_frame, p_U_out_frame, p_abutting = isempty(p_gaps) ? (0, 0, 0, 0, 0, 0, 0, 0) : count_overlaps(p_gaps)
                n_U_overlap_nr, n_U_overlap_length_sum, n_U_max_len, n_U_coupled_nr, n_U_deep_nr, n_U_in_frame, n_U_out_frame, n_abutting = isempty(n_gaps) ? (0, 0, 0, 0, 0, 0, 0, 0) : count_overlaps(n_gaps)
                
                p_gap_length_sum = isempty(p_gaps) ? 0 : sum(filter(x -> x > 0, p_gaps))
                n_gap_length_sum = isempty(n_gaps) ? 0 : sum(filter(x -> x > 0, n_gaps))
                
                if !isempty(p_genes) && !isempty(n_genes)
                    p_intervals = GenomicFeatures.Interval.(Ref(contig_id), p_starts, p_ends, Ref(STRAND_POS), 1:p_gene_nr)
                    n_intervals = GenomicFeatures.Interval.(Ref(contig_id), n_starts, n_ends, Ref(STRAND_NEG), 1:n_gene_nr)
                    C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum, C_max_len, D_max_len, C_deep_nr, D_deep_nr = bidirectional_overlaps(
                        IntervalCollection(n_intervals, true), 
                        IntervalCollection(p_intervals, true)
                    )
                else
                    C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum, C_max_len, D_max_len, C_deep_nr, D_deep_nr = 0, 0, 0, 0, 0, 0, 0, 0
                end
                
                p_gap_stats = analyze_intergenic_distances(p_gaps)
                n_gap_stats = analyze_intergenic_distances(n_gaps)
                operon_stats = analyze_operons(all_genes_sorted)
                strand_asymmetry = isempty(genes) ? 0.5 : p_gene_nr / num_genes
                local_arrangements = analyze_local_arrangements(all_genes_sorted)
                density_gradient_std = calculate_density_gradient(all_genes_sorted, contig_size)
                
                coding_density_pct = calculate_coding_density(all_genes_sorted, contig_size)
                nested_genes_nr = count_nested_genes(all_genes_sorted)

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
                    
                    absolute_gap_mean=absolute_gap_stats.mean, 
                    absolute_gap_median=absolute_gap_stats.median, 
                    absolute_gap_std=absolute_gap_stats.std,

                    operon_nr=operon_stats.operon_nr, operonicity_score=operon_stats.operonicity_score,
                    mean_operon_size=operon_stats.mean_operon_size, strand_asymmetry=strand_asymmetry,
                    divergent_pairs_nr=local_arrangements.divergent_pairs_nr,
                    convergent_pairs_nr=local_arrangements.convergent_pairs_nr,
                    gene_density_gradient_std=density_gradient_std,
                    
                    U_coupled_nr = p_U_coupled_nr + n_U_coupled_nr,
                    U_deep_nr = p_U_deep_nr + n_U_deep_nr,
                    U_max_len = max(p_U_max_len, n_U_max_len),
                    C_deep_nr = C_deep_nr,
                    D_deep_nr = D_deep_nr,
                    C_max_len = C_max_len,
                    D_max_len = D_max_len,
                    
                    U_in_frame_nr = p_U_in_frame + n_U_in_frame,
                    U_out_of_frame_nr = p_U_out_frame + n_U_out_frame,
                    abutting_genes_nr = p_abutting + n_abutting,
                    strand_switch_rate = strand_switch_rate,
                    mean_gene_length = mean_gene_length,
                    std_gene_length = std_gene_length,
                    coding_density_pct = coding_density_pct,
                    nested_genes_nr = nested_genes_nr,
                    
                    mean_overlap_chain_size = overlap_chains.mean_chain_size,
                    max_overlap_chain_size = overlap_chains.max_chain_size,
                    total_overlap_chains_nr = overlap_chains.total_chains,
                    intra_operon_gap_mean = accordion_gaps.intra_gap_mean,
                    inter_operon_gap_mean = accordion_gaps.inter_gap_mean,
                    max_overlaps_per_gene = entanglement_hubs_max,
                    
                    internal_overlaps_nr = p_U_overlap_nr + n_U_overlap_nr,
                    boundary_collisions_nr = C_overlap_nr + D_overlap_nr
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