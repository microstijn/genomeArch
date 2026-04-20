using CSV
using DataFrames
using BioSequences
using FASTX
using GenomicFeatures
using CairoMakie
using GFF3
println("Initializing the Master Hexamer Training & GHMM Pipeline...")

# ==============================================================================
# 1. HEXAMER TRAINING ENGINE (From your script)
# ==============================================================================
function initialize_hexamer_dict()
    bases = [dna"A", dna"C", dna"G", dna"T"]
    dict = Dict{LongDNA{4}, Float64}()
    for b1 in bases, b2 in bases, b3 in bases, b4 in bases, b5 in bases, b6 in bases
        dict[b1 * b2 * b3 * b4 * b5 * b6] = 1.0 
    end
    return dict
end

function parse_gff_for_cds(gff_path::String)
    cds_coords = Tuple{Int, Int, Char}[]
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) == "CDS"
            start_pos = GFF3.seqstart(record)
            end_pos = GFF3.seqend(record)
            strand_val = GFF3.strand(record)
            strand_char = strand_val == GenomicFeatures.Strand('+') ? '+' : '-'
            push!(cds_coords, (start_pos, end_pos, strand_char))
        end
    end
    close(reader)
    return cds_coords
end

function train_hexamer_matrices(fasta_path::String, gff_path::String)
    println("  -> Loading Genome for Training...")
    reader = FASTA.Reader(open(fasta_path, "r"))
    genome = LongDNA{4}(FASTA.sequence(first(reader)))
    close(reader)
    
    println("  -> Parsing GFF & Building Frequencies...")
    cds_coords = parse_gff_for_cds(gff_path)
    
    counts_phase0 = initialize_hexamer_dict(); counts_phase1 = initialize_hexamer_dict()
    counts_phase2 = initialize_hexamer_dict(); counts_background = initialize_hexamer_dict()
    total_p0 = 4096.0; total_p1 = 4096.0; total_p2 = 4096.0; total_bg = 4096.0
    
    is_coding = falses(length(genome)) 
    
    for (start_pos, end_pos, strand) in cds_coords
        # --- THE FIX: Clamp the coordinates to the physical array size ---
        safe_start = max(1, start_pos)
        safe_end = min(length(genome), end_pos)
        
        # If a plasmid annotation snuck in and is completely out of bounds, skip it
        if safe_start >= safe_end continue end 
        
        is_coding[safe_start:safe_end] .= true
        seq = genome[safe_start:safe_end]
        if strand == '-' seq = reverse_complement(seq) end
        
        for i in 1:(length(seq)-5)
            hexamer = seq[i:i+5]
            if hasambiguity(hexamer) continue end
            phase = (i - 1) % 3
            if phase == 0     counts_phase0[hexamer] += 1.0; total_p0 += 1.0
            elseif phase == 1 counts_phase1[hexamer] += 1.0; total_p1 += 1.0
            else              counts_phase2[hexamer] += 1.0; total_p2 += 1.0 end
        end
    end
    
    for i in 1:(length(genome)-5)
        if !any(is_coding[i:i+5])
            hexamer = genome[i:i+5]
            if hasambiguity(hexamer) continue end
            counts_background[hexamer] += 1.0; total_bg += 1.0
            counts_background[reverse_complement(hexamer)] += 1.0; total_bg += 1.0
        end
    end
    
    matrix_p0 = Dict{LongDNA{4}, Float64}()
    matrix_p1 = Dict{LongDNA{4}, Float64}()
    matrix_p2 = Dict{LongDNA{4}, Float64}()
    
    for (hex, bg_count) in counts_background
        prob_bg = bg_count / total_bg
        matrix_p0[hex] = log2((counts_phase0[hex] / total_p0) / prob_bg)
        matrix_p1[hex] = log2((counts_phase1[hex] / total_p1) / prob_bg)
        matrix_p2[hex] = log2((counts_phase2[hex] / total_p2) / prob_bg)
    end
    
    println("  -> GHMM Matrices Successfully Trained.")
    return (matrix_p0, matrix_p1, matrix_p2)
end

# ==============================================================================
# 2. GHMM SCORING & EXTRACTION
# ==============================================================================
function window_ghmm_score(win_seq::LongDNA{4}, matrices, absolute_start::Int)
    score = 0.0
    hexamer_count = length(win_seq) - 5
    if hexamer_count <= 0 return 0.0 end
    
    for j in 1:hexamer_count
        hexamer = win_seq[j:j+5]
        if hasambiguity(hexamer) continue end
        
        absolute_pos = absolute_start + j - 1
        phase = (absolute_pos - 1) % 3
        
        if phase == 0     score += get(matrices[1], hexamer, -0.2)
        elseif phase == 1 score += get(matrices[2], hexamer, -0.4)
        else              score += get(matrices[3], hexamer, -0.4) end
    end
    return score / hexamer_count
end

function slide_dual_ghmm(seq_fwd::LongDNA{4}, matrices, abs_start::Int; win_nt=90)
    positions = Int[]; scores_fwd = Float64[]; scores_rev = Float64[]
    seq_rev = reverse_complement(seq_fwd)
    seq_len = length(seq_fwd)
    
    for i in 1:(seq_len - win_nt + 1)
        # Forward window
        subseq_fwd = seq_fwd[i:i+win_nt-1]
        score_f = window_ghmm_score(subseq_fwd, matrices, abs_start + i - 1)
        
        # Reverse window
        rev_start = seq_len - (i + win_nt - 1) + 1
        subseq_rev = seq_rev[rev_start:rev_start+win_nt-1]
        score_r = window_ghmm_score(subseq_rev, matrices, rev_start)
        
        push!(positions, i + (win_nt ÷ 2))
        push!(scores_fwd, score_f)
        push!(scores_rev, score_r)
    end
    return positions, scores_fwd, scores_rev
end

function get_locus_data(gff_path, fasta_path, gene_A_id, gene_B_id, pad=300)
    start_A, end_A, start_B, end_B = -1, -1, -1, -1
    
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) ∈ ("CDS", "gene")
            attr_str = string(GFF3.attributes(record))
            if occursin(gene_A_id, attr_str)
                start_A, end_A = GFF3.seqstart(record), GFF3.seqend(record)
            elseif occursin(gene_B_id, attr_str)
                start_B, end_B = GFF3.seqstart(record), GFF3.seqend(record)
            end
        end
    end
    close(reader)
    
    if start_A == -1 || start_B == -1 return nothing, -1, -1, -1, -1 end
    
    locus_start = max(1, min(start_A, start_B) - pad)
    overlap_start = max(start_A, start_B)
    overlap_end = min(end_A, end_B)
    
    seq = dna""
    f_reader = FASTA.Reader(open(fasta_path, "r"))
    for record in f_reader
        full_seq = FASTA.sequence(record)
        locus_end = min(length(full_seq), max(end_A, end_B) + pad)
        seq = full_seq[locus_start:locus_end]
        break
    end
    close(f_reader)
    
    rel_overlap_start = overlap_start - locus_start + 1
    rel_overlap_end = overlap_end - locus_start + 1
    
    return seq, rel_overlap_start, rel_overlap_end, locus_start, (locus_start + length(seq) - 1)
end

# ==============================================================================
# 3. PLOTTING ENGINE
# ==============================================================================
function plot_dual(positions, s_fwd, s_rev, overlap_start, overlap_end, title_str, out_path)
    set_theme!(theme_minimal(), font="Helvetica")
    fig = Figure(size = (900, 500))
    ax = Axis(fig[1, 1], title = title_str, xlabel = "Relative Locus Position (bp)", ylabel = "GHMM Log-Odds Score")
    
    vspan!(ax, overlap_start, overlap_end, color = (:goldenrod, 0.2), label = "Physical Overlap Zone")
    hlines!(ax, [0.0], color = :black, linestyle = :dash, label = "Coding Threshold")
    
    lines!(ax, positions, s_fwd, color = :midnightblue, linewidth = 2.5, label = "Forward Frame Score")
    lines!(ax, positions, s_rev, color = :crimson, linewidth = 2.5, label = "Reverse Frame Score")
    
    axislegend(ax, position = :rt, framevisible = false)
    save(out_path, fig)
end

# ==============================================================================
# 4. MAIN EXECUTION LOOP
# ==============================================================================
# Load the target dataset
survivors_df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)

# Group the overlaps by genome to avoid re-training hexamers
grouped_genomes = groupby(survivors_df, :genome)

out_dir = raw"D:\pipeline_output\trained_ghmm_plots"
isdir(out_dir) || mkdir(out_dir)
base_ncbi_dir = raw"D:\ncbi_downloads\target_genomes\survivor_data\ncbi_dataset\data"

println("Starting processing of $(length(grouped_genomes)) unique genomes...")

for sub_df in grouped_genomes
    genome_acc = sub_df.genome[1]
    println("\n========================================")
    println("Processing Genome: $genome_acc")
    
    genome_dir = joinpath(base_ncbi_dir, genome_acc)
    if !isdir(genome_dir)
        println("⚠️  Skipping: Directory not found at $genome_dir")
        continue
    end
    
    fasta_files = filter(f -> endswith(f, ".fna"), readdir(genome_dir))
    gff_files   = filter(f -> endswith(f, ".gff"), readdir(genome_dir))
    
    if isempty(fasta_files) || isempty(gff_files)
        println("⚠️  Skipping: Missing .fna or .gff in $genome_dir")
        continue
    end
    
    fasta_path = joinpath(genome_dir, fasta_files[1])
    gff_path   = joinpath(genome_dir, gff_files[1])
    
    # 1. Train the Hexamers for this specific genome
    println("Training empirical GHMM matrices...")
    matrices = train_hexamer_matrices(fasta_path, gff_path)
    
    # 2. Iterate through all overlaps in this genome
    println("Evaluating $(nrow(sub_df)) overlapping loci...")
    for row in eachrow(sub_df)
        clean_gene_A = replace(row.gene_A, "gene-" => "")
        clean_gene_B = replace(row.gene_B, "gene-" => "")
        
        seq, rel_start, rel_end, abs_start, abs_end = get_locus_data(gff_path, fasta_path, clean_gene_A, clean_gene_B)
        
        if seq === nothing
            println("  -> ⚠️  Skipped: Coords not found for $clean_gene_A / $clean_gene_B")
            continue
        end
        
        # Calculate real GHMM scores based on our trained matrices
        # We explicitly wrap `seq` in LongDNA{4}() to guarantee the type matches
        pos, s_fwd, s_rev = slide_dual_ghmm(LongDNA{4}(seq), matrices, abs_start, win_nt=90)
        
        # Plot
        safe_name = replace("$(genome_acc)_$(clean_gene_A)_vs_$(clean_gene_B)", r"[^a-zA-Z0-9_]" => "_")
        plot_name = joinpath(out_dir, "$(safe_name).pdf")
        
        plot_dual(pos, s_fwd, s_rev, rel_start, rel_end, 
                  "Locus: $genome_acc | Overlap: $(row.overlap_length)bp", plot_name)
        
        println("  -> ✅ Plotted: $(clean_gene_A) vs $(clean_gene_B)")
    end
end

println("\nComplete! All empirical GHMM trajectories saved to '$out_dir'.")