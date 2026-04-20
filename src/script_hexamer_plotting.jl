using CSV
using DataFrames
using BioSequences
using FASTX
using GenomicFeatures
using CairoMakie
using JSON3

println("🚀 Initializing the Ultimate Offline Domain Plotting Pipeline...")

# ==============================================================================
# 1. STRUCTURAL DOMAIN CACHE LOADING
# ==============================================================================
struct PlotDomain
    name::String
    start_aa::Int
    end_aa::Int
    color::Symbol
end

println("Loading pre-computed structural domain cache...")
raw_cache = read(raw"D:\pipeline_output\domain_cache.json", String)
domain_cache = Dict{String, Vector{PlotDomain}}()

for (k, v) in JSON3.read(raw_cache)
    # Convert the saved string color back into a Julia Symbol for CairoMakie
    domain_cache[String(k)] = [PlotDomain(d.name, d.start_aa, d.end_aa, Symbol(d.color)) for d in v]
end

function map_domain_to_dna(domain::PlotDomain, start_codon::Int, strand::Char)
    if strand == '+'
        dna_start = start_codon + (domain.start_aa - 1) * 3
        dna_end   = start_codon + (domain.end_aa * 3) - 1
        return (dna_start, dna_end)
    else
        dna_start = start_codon - (domain.end_aa * 3) + 1
        dna_end   = start_codon - (domain.start_aa - 1) * 3
        return (dna_start, dna_end)
    end
end

# ==============================================================================
# 2. 5TH-ORDER MARKOV TRAINING ENGINE
# ==============================================================================
function initialize_hexamer_dict()
    bases = [dna"A", dna"C", dna"G", dna"T"]
    dict = Dict{LongDNA{4}, Float64}()
    for b1 in bases, b2 in bases, b3 in bases, b4 in bases, b5 in bases, b6 in bases
        dict[b1 * b2 * b3 * b4 * b5 * b6] = 1.0 
    end
    return dict
end

function initialize_pentamer_dict()
    bases = [dna"A", dna"C", dna"G", dna"T"]
    dict = Dict{LongDNA{4}, Float64}()
    for b1 in bases, b2 in bases, b3 in bases, b4 in bases, b5 in bases
        dict[b1 * b2 * b3 * b4 * b5] = 4.0 
    end
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

# ==============================================================================
# 3. LOCUS EXTRACTION & GHMM SLIDER
# ==============================================================================
function get_locus_data(gff_path, fasta_path, gene_A_id, gene_B_id, pad=300)
    start_A, end_A, start_B, end_B = -1, -1, -1, -1
    strand_A, strand_B = '+', '+' 
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) ∈ ("CDS", "gene")
            attr_str = string(GFF3.attributes(record))
            if occursin(gene_A_id, attr_str)
                start_A, end_A = GFF3.seqstart(record), GFF3.seqend(record)
                strand_A = GFF3.strand(record) == GenomicFeatures.Strand('+') ? '+' : '-'
            elseif occursin(gene_B_id, attr_str)
                start_B, end_B = GFF3.seqstart(record), GFF3.seqend(record)
                strand_B = GFF3.strand(record) == GenomicFeatures.Strand('+') ? '+' : '-'
            end
        end
    end
    close(reader)
    if start_A == -1 || start_B == -1 return nothing end
    
    locus_start = max(1, min(start_A, start_B) - pad)
    
    seq = dna""
    f_reader = FASTA.Reader(open(fasta_path, "r"))
    for record in f_reader
        full_seq = LongDNA{4}(FASTA.sequence(record))
        seq = full_seq[locus_start:min(length(full_seq), max(end_A, end_B) + pad)]
        break
    end
    close(f_reader)
    
    return (seq=seq, rel_start=max(start_A, start_B) - locus_start + 1, rel_end=min(end_A, end_B) - locus_start + 1, 
            abs_start=locus_start, start_A=start_A, end_A=end_A, strand_A=strand_A, 
            start_B=start_B, end_B=end_B, strand_B=strand_B)
end

function slide_frame_locked_ghmm(seq_fwd::LongDNA{4}, matrices, abs_start::Int, gene_A_start::Int, gene_B_start::Int; win_nt=90)
    positions = Int[]; scores_fwd = Float64[]; scores_rev = Float64[]
    seq_rev = reverse_complement(seq_fwd)
    seq_len = length(seq_fwd)
    
    for i in 1:(seq_len - win_nt + 1)
        slice_genomic_start = abs_start + i - 1
        
        score_f = 0.0; count_f = 0
        for j in 1:(win_nt-5)
            hex = seq_fwd[i+j-1:i+j+4]
            if !hasambiguity(hex)
                phase = mod((slice_genomic_start + j - 1) - gene_A_start, 3)
                score_f += get(matrices[phase+1], hex, -0.2); count_f += 1
            end
        end
        
        score_r = 0.0; count_r = 0
        rev_start = seq_len - (i + win_nt - 1) + 1
        for j in 1:(win_nt-5)
            hex = seq_rev[rev_start+j-1:rev_start+j+4]
            if !hasambiguity(hex)
                phase = mod(gene_B_start - ((slice_genomic_start + win_nt - 1) - j + 1), 3)
                score_r += get(matrices[phase+1], hex, -0.4); count_r += 1
            end
        end
        
        push!(positions, i + (win_nt ÷ 2))
        push!(scores_fwd, count_f > 0 ? score_f/count_f : 0.0)
        push!(scores_rev, count_r > 0 ? score_r/count_r : 0.0)
    end
    return positions, scores_fwd, scores_rev
end

# ==============================================================================
# 4. DOMAIN-AWARE CAIROMAKIE PLOTTING
# ==============================================================================
function plot_biological_with_domains(positions, s_fwd, s_rev, locus, domains_A, domains_B, title_str, row, out_path)
    dir_A = row.strand_A == "+" ? "Left → Right" : "Right ← Left"
    dir_B = row.strand_B == "+" ? "Left → Right" : "Right ← Left"
    sub_str = "Type: $(row.overlap_type) | Host ($(row.strand_A)): $dir_A | Antisense ($(row.strand_B)): $dir_B"

    set_theme!(theme_minimal(), font="Helvetica")
    
    # Increase height to comfortably fit two panels
    fig = Figure(size = (1000, 800))
    
    # --- TOP PANEL (Log-Odds Scores) ---
    ax_top = Axis(fig[1, 1], 
        title = title_str, 
        subtitle = sub_str, 
        ylabel = "5th-Order Markov Score"
    )
    # Hide X-axis labels on the top plot to prevent clutter
    hidexdecorations!(ax_top, grid = false) 
    
    # --- BOTTOM PANEL (Structural Domains) ---
    ax_bot = Axis(fig[2, 1], 
        xlabel = "Relative Locus Position (bp)", 
        ylabel = "Domains",
        # Create custom clean Y-ticks for the two gene tracks
        yticks = ([-1.0, 1.0], ["Gene B (-)", "Gene A (+)"])
    )
    
    # Link the X-axes so they share the same physical coordinates
    linkxaxes!(ax_top, ax_bot)
    
    # Make the top panel take up 70% of the height, bottom takes 30%
    rowsize!(fig.layout, 1, Relative(0.7))

    # --- 1. NaN MASKING (Math) ---
    masked_fwd = [(locus.abs_start + p - 1 >= locus.start_A && locus.abs_start + p - 1 <= locus.end_A) ? s_fwd[i] : NaN for (i,p) in enumerate(positions)]
    masked_rev = [(locus.abs_start + p - 1 >= locus.start_B && locus.abs_start + p - 1 <= locus.end_B) ? s_rev[i] : NaN for (i,p) in enumerate(positions)]

    # --- 2. SHARED VISUALS (X-Limits & Overlap Shading) ---
    crop_start = min(locus.start_A, locus.start_B) - locus.abs_start + 1
    crop_end   = max(locus.end_A, locus.end_B) - locus.abs_start + 1
    xlims!(ax_top, crop_start, crop_end)
    
    # Shade the physical overlap on BOTH panels
    # --- SHARED VISUALS (Overlap vs Gap Shading) ---
    if locus.rel_start <= locus.rel_end
        # TRUE OVERLAP (Goldenrod)
        vspan!(ax_top, locus.rel_start, locus.rel_end, color = (:goldenrod, 0.15), label = "Physical Overlap")
        vspan!(ax_bot, locus.rel_start, locus.rel_end, color = (:goldenrod, 0.15))
    else
        # INTERGENIC GAP (Gray / Red)
        # We swap the start and end so it draws left-to-right correctly
        vspan!(ax_top, locus.rel_end, locus.rel_start, color = (:lightcoral, 0.15), label = "Intergenic Gap")
        vspan!(ax_bot, locus.rel_end, locus.rel_start, color = (:lightcoral, 0.15))
    end

    hlines!(ax_top, [0.0], color = :black, linestyle = :solid, linewidth = 1.0)
    hlines!(ax_bot, [0.0], color = :gray80, linestyle = :dash, linewidth = 1.0)

    # --- 3. PLOT LOG-ODDS LINES (Top Panel) ---
    lines!(ax_top, positions, masked_fwd, color = :midnightblue, linewidth = 2.5, label = "Gene A Score")
    lines!(ax_top, positions, masked_rev, color = :crimson, linewidth = 2.5, label = "Gene B Score")
    axislegend(ax_top, position = :rt, framevisible = true)

    # --- 4. DOMAIN TRACKS (Bottom Panel) ---
    # Gene A domains sit high at Y = 1.0
    true_start_A = row.strand_A == "+" ? locus.start_A : locus.end_A
    for dom in domains_A
        d_start, d_end = map_domain_to_dna(dom, true_start_A, row.strand_A[1])
        r_start, r_end = d_start - locus.abs_start + 1, d_end - locus.abs_start + 1
        
        # Force color based on Gene A
        lines!(ax_bot, [r_start, r_end], [1.0, 1.0], color=:midnightblue, linewidth=20)
        text!(ax_bot, (r_start + r_end)/2, 1.0, text=dom.name, color=:white, align=(:center, :center), fontsize=11, font="Helvetica Bold")
    end
    
    # Gene B domains sit low at Y = -1.0
    true_start_B = row.strand_B == "+" ? locus.start_B : locus.end_B
    for dom in domains_B
        d_start, d_end = map_domain_to_dna(dom, true_start_B, row.strand_B[1])
        r_start, r_end = d_start - locus.abs_start + 1, d_end - locus.abs_start + 1
        
        # Force color based on Gene B
        lines!(ax_bot, [r_start, r_end], [-1.0, -1.0], color=:crimson, linewidth=20)
        text!(ax_bot, (r_start + r_end)/2, -1.0, text=dom.name, color=:white, align=(:center, :center), fontsize=11, font="Helvetica Bold")
    end
    
    ylims!(ax_bot, -2.0, 2.0)

    # Save as PNG
    save(out_path, fig, px_per_unit = 2) # High-resolution PNG scaling
end
# ==============================================================================
# 5. MASTER EXECUTION LOOP
# ==============================================================================
survivors_df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)

grouped_genomes = groupby(survivors_df, :genome)
base_ncbi_dir = raw"D:\ncbi_downloads\target_genomes\survivor_data\ncbi_dataset\data"

out_dir = raw"D:\pipeline_output\structural_domain_plots"
isdir(out_dir) || mkdir(out_dir)

println("\nStarting Local Genome Processing Loop...")

unique(survivors_df.genome) |> println

for sub_df in grouped_genomes
    genome_acc = sub_df.genome[1]
    genome_dir = joinpath(base_ncbi_dir, genome_acc)
    
    if !isdir(genome_dir) continue end
    fasta_files = filter(f -> endswith(f, ".fna"), readdir(genome_dir))
    gff_files   = filter(f -> endswith(f, ".gff"), readdir(genome_dir))
    if isempty(fasta_files) || isempty(gff_files) continue end
    
    println("\n-> Training 5th-Order Matrices for $genome_acc...")
    matrices = train_5th_order_matrices(joinpath(genome_dir, fasta_files[1]), joinpath(genome_dir, gff_files[1]))
    
    for row in eachrow(sub_df)
        clean_gene_A = replace(row.gene_A, "gene-" => "")
        clean_gene_B = replace(row.gene_B, "gene-" => "")
        
        locus = get_locus_data(joinpath(genome_dir, gff_files[1]), joinpath(genome_dir, fasta_files[1]), clean_gene_A, clean_gene_B)
        if locus === nothing continue end
        
        # Fetch the parsed domains directly from the JSON dictionary!
        domains_A = get(domain_cache, clean_gene_A, PlotDomain[])
        domains_B = get(domain_cache, clean_gene_B, PlotDomain[])
        
        true_start_A = locus.strand_A == '+' ? locus.start_A : locus.end_A
        true_start_B = locus.strand_B == '+' ? locus.start_B : locus.end_B
        pos, s_fwd, s_rev = slide_frame_locked_ghmm(LongDNA{4}(locus.seq), matrices, locus.abs_start, true_start_A, true_start_B, win_nt=90)
        
        plot_name = joinpath(out_dir, "$(genome_acc)_$(clean_gene_A)_vs_$(clean_gene_B).png")
        
        # Dynamically build a title
        domain_status = "Domains: A($(length(domains_A))) | B($(length(domains_B)))"
        title = "Locus: $genome_acc | Overlap: $(row.overlap_length)bp | $domain_status\n$clean_gene_A vs $clean_gene_B"
        
        plot_biological_with_domains(pos, s_fwd, s_rev, locus, domains_A, domains_B, title, row, plot_name)
        println("   ✅ Generated: $plot_name")
    end
end

println("\n🎉 Done! All 97 overlap visualizations are complete. Check the '$out_dir' directory.")