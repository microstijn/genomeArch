using CSV
using DataFrames
using BioSequences
using FASTX
using GenomicFeatures
using HypothesisTests

println("🧮 Initializing Combined Paired Markov Statistical Analysis...")

# ==============================================================================
# 1. MARKOV SCORING ENGINE
# ==============================================================================
function initialize_hexamer_dict()
    bases = [dna"A", dna"C", dna"G", dna"T"]
    dict = Dict{LongDNA{4}, Float64}()
    for b1 in bases, b2 in bases, b3 in bases, b4 in bases, b5 in bases, b6 in bases dict[b1*b2*b3*b4*b5*b6] = 1.0 end
    return dict
end
function initialize_pentamer_dict()
    bases = [dna"A", dna"C", dna"G", dna"T"]
    dict = Dict{LongDNA{4}, Float64}()
    for b1 in bases, b2 in bases, b3 in bases, b4 in bases, b5 in bases dict[b1*b2*b3*b4*b5] = 4.0 end
    return dict
end

function train_5th_order_matrices(fasta_path::String, gff_path::String)
    reader = FASTA.Reader(open(fasta_path, "r"))
    genome = LongDNA{4}(FASTA.sequence(first(reader))) 
    close(reader)
    
    cds_coords = Tuple{Int, Int, Char}[]
    gff_reader = GFF3.Reader(open(gff_path, "r"))
    for record in gff_reader
        if GFF3.featuretype(record) == "CDS"
            strand = GFF3.strand(record) == GenomicFeatures.Strand('+') ? '+' : '-'
            push!(cds_coords, (GFF3.seqstart(record), GFF3.seqend(record), strand))
        end
    end
    close(gff_reader)
    
    counts_hex_p0 = initialize_hexamer_dict(); counts_hex_p1 = initialize_hexamer_dict()
    counts_hex_p2 = initialize_hexamer_dict(); counts_hex_bg = initialize_hexamer_dict()
    counts_pent_p0 = initialize_pentamer_dict(); counts_pent_p1 = initialize_pentamer_dict()
    counts_pent_p2 = initialize_pentamer_dict(); counts_pent_bg = initialize_pentamer_dict()
    
    is_coding = falses(length(genome)) 
    for (start_pos, end_pos, strand) in cds_coords
        safe_start, safe_end = max(1, start_pos), min(length(genome), end_pos)
        if safe_start >= safe_end continue end 
        is_coding[safe_start:safe_end] .= true
        seq = strand == '-' ? reverse_complement(genome[safe_start:safe_end]) : genome[safe_start:safe_end]
        for i in 1:(length(seq)-5)
            hexamer = seq[i:i+5]
            if hasambiguity(hexamer) continue end
            prefix, phase = hexamer[1:5], (i - 1) % 3
            if phase == 0     counts_hex_p0[hexamer] += 1.0; counts_pent_p0[prefix] += 1.0
            elseif phase == 1 counts_hex_p1[hexamer] += 1.0; counts_pent_p1[prefix] += 1.0
            else              counts_hex_p2[hexamer] += 1.0; counts_pent_p2[prefix] += 1.0 end
        end
    end
    for i in 1:(length(genome)-5)
        if !any(is_coding[i:i+5])
            hexamer = genome[i:i+5]
            if hasambiguity(hexamer) continue end
            prefix = hexamer[1:5]
            counts_hex_bg[hexamer] += 1.0; counts_pent_bg[prefix] += 1.0
            rc_hexamer = reverse_complement(hexamer)
            counts_hex_bg[rc_hexamer] += 1.0; counts_pent_bg[rc_hexamer[1:5]] += 1.0
        end
    end
    
    matrix_p0 = Dict{LongDNA{4}, Float64}(); matrix_p1 = Dict{LongDNA{4}, Float64}(); matrix_p2 = Dict{LongDNA{4}, Float64}()
    for (hex, bg_hex_count) in counts_hex_bg
        prefix = hex[1:5]
        prob_bg = bg_hex_count / counts_pent_bg[prefix]
        matrix_p0[hex] = log2((counts_hex_p0[hex] / counts_pent_p0[prefix]) / prob_bg)
        matrix_p1[hex] = log2((counts_hex_p1[hex] / counts_pent_p1[prefix]) / prob_bg)
        matrix_p2[hex] = log2((counts_hex_p2[hex] / counts_pent_p2[prefix]) / prob_bg)
    end
    return (matrix_p0, matrix_p1, matrix_p2)
end

function score_region(seq::LongDNA{4}, matrices)
    if length(seq) < 6 return missing end
    score_sum = 0.0
    count = 0
    for i in 1:(length(seq)-5)
        hex = seq[i:i+5]
        if !hasambiguity(hex)
            phase = (i - 1) % 3
            score_sum += get(matrices[phase+1], hex, 0.0)
            count += 1
        end
    end
    return count > 0 ? score_sum / count : missing
end

# ==============================================================================
# 2. SEQUENCE EXTRACTION & DATA AGGREGATION
# ==============================================================================
survivors_df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)
base_ncbi_dir = raw"D:\ncbi_downloads\target_genomes\survivor_data\ncbi_dataset\data"

grouped_genomes = groupby(survivors_df, :genome)


# Removed "Role" - everything is grouped together
results = DataFrame(
    Genome = String[], Gene = String[], Overlap_Type = String[],
    Mean_Outside = Float64[], Mean_Inside = Float64[], Delta = Float64[], Category = String[]
)

for sub_df in grouped_genomes
    genome_acc = sub_df.genome[1]
    genome_dir = joinpath(base_ncbi_dir, genome_acc)
    fasta_files = filter(f -> endswith(f, ".fna"), readdir(genome_dir))
    gff_files   = filter(f -> endswith(f, ".gff") || endswith(f, ".gff3"), readdir(genome_dir))
    if isempty(fasta_files) || isempty(gff_files) continue end
    
    matrices = train_5th_order_matrices(joinpath(genome_dir, fasta_files[1]), joinpath(genome_dir, gff_files[1]))
    
    reader = FASTA.Reader(open(joinpath(genome_dir, fasta_files[1]), "r"))
    genome_seq = LongDNA{4}(FASTA.sequence(first(reader)))
    close(reader)
    
    gff_path = joinpath(genome_dir, gff_files[1])
    
    for row in eachrow(sub_df)
        clean_gene_A = replace(row.gene_A, "gene-" => "")
        clean_gene_B = replace(row.gene_B, "gene-" => "")
        
        function get_coords(gene_id)
            r = GFF3.Reader(open(gff_path, "r"))
            for rec in r
                if GFF3.featuretype(rec) == "CDS" && occursin(gene_id, string(GFF3.attributes(rec)))
                    start_p, end_p = GFF3.seqstart(rec), GFF3.seqend(rec)
                    strand = GFF3.strand(rec) == GenomicFeatures.Strand('+') ? '+' : '-'
                    close(r); return (start_p, end_p, strand)
                end
            end
            close(r); return nothing
        end
        
        coords_A = get_coords(clean_gene_A)
        coords_B = get_coords(clean_gene_B)
        if coords_A === nothing || coords_B === nothing continue end
        
        overlap_start = max(coords_A[1], coords_B[1])
        overlap_end = min(coords_A[2], coords_B[2])
        if overlap_start > overlap_end continue end # Skip false gaps
        
        # Slices, orients, and scores a gene
        function process_gene(coords, gene_id)
            g_start, g_end, strand = coords
            
            full_seq = strand == '+' ? genome_seq[g_start:g_end] : reverse_complement(genome_seq[g_start:g_end])
            
            rel_over_start = strand == '+' ? (overlap_start - g_start + 1) : (g_end - overlap_end + 1)
            rel_over_end   = strand == '+' ? (overlap_end - g_start + 1) : (g_end - overlap_start + 1)
            
            seq_inside = full_seq[rel_over_start:rel_over_end]
            seq_outside = full_seq[1:rel_over_start-1] * full_seq[rel_over_end+1:end]
            
            score_in = score_region(seq_inside, matrices)
            score_out = score_region(seq_outside, matrices)
            
            if !ismissing(score_in) && !ismissing(score_out)
                delta = score_in - score_out
                cat = delta > 0 ? "Higher" : (delta < 0 ? "Lower" : "Same")
                push!(results, (genome_acc, gene_id, row.overlap_type, score_out, score_in, delta, cat))
            end
        end
        
        process_gene(coords_A, clean_gene_A)
        process_gene(coords_B, clean_gene_B)
    end
end

CSV.write("markov_combined_statistics.csv", results)
println("✅ Saved raw paired data to 'markov_combined_statistics.csv'")

# ==============================================================================
# 3. STATISTICAL TESTS
# ==============================================================================
higher_count = sum(results.Category .== "Higher")
lower_count = sum(results.Category .== "Lower")
total = higher_count + lower_count # Excludes the extremely rare "Same" ties

println("\n" * "="^50)
println("🔬 COMBINED STATISTICAL RESULTS: OVERLAP VS NON-OVERLAP")
println("="^50)
println("Total Genes Evaluated: $(nrow(results))")
println("Overlap Score HIGHER:  $higher_count")
println("Overlap Score LOWER:   $lower_count")

# 1. The Categorical Sign Test
bt = BinomialTest(higher_count, total, 0.5)
println("\n1. Categorical Sign Test (Null: 50% Higher / 50% Lower)")
println("   p-value = $(round(pvalue(bt), digits=6))")

# 2. Paired Wilcoxon Signed-Rank Test
wt = SignedRankTest(results.Mean_Inside, results.Mean_Outside)
println("\n2. Paired Wilcoxon Signed-Rank Test (Null: Median Difference is 0)")
println("   p-value = $(round(pvalue(wt), digits=6))")
println("="^50)