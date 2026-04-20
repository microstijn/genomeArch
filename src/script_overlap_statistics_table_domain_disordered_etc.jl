using CSV
using DataFrames
using JSON3
using GenomicFeatures

println("📊 Initializing Structural Overlap Statistical Analysis...")

# ==============================================================================
# 1. SETUP & CACHE LOADING
# ==============================================================================
struct PlotDomain
    name::String
    start_aa::Int
    end_aa::Int
    color::Symbol
end

raw_cache = JSON3.read(read("domain_cache.json", String))
domain_cache = Dict{String, Vector{PlotDomain}}()
for (k, v) in raw_cache
    domain_cache[String(k)] = [PlotDomain(d.name, d.start_aa, d.end_aa, Symbol(d.color)) for d in v]
end

function map_domain_to_dna(domain::PlotDomain, start_codon::Int, strand::Char)
    if strand == '+'
        return (start_codon + (domain.start_aa - 1) * 3, start_codon + (domain.end_aa * 3) - 1)
    else
        return (start_codon - (domain.end_aa * 3) + 1, start_codon - (domain.start_aa - 1) * 3)
    end
end

function get_gene_coords(gff_path, gene_id)
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) ∈ ("CDS", "gene")
            if occursin(gene_id, string(GFF3.attributes(record)))
                start_pos = GFF3.seqstart(record)
                end_pos = GFF3.seqend(record)
                strand = GFF3.strand(record) == GenomicFeatures.Strand('+') ? '+' : '-'
                close(reader)
                return (start_pos, end_pos, strand)
            end
        end
    end
    close(reader)
    return nothing
end

# ==============================================================================
# 2. OVERLAP CLASSIFIER
# ==============================================================================
function classify_overlap_zone(domains::Vector{PlotDomain}, start_codon, strand, overlap_start, overlap_end)
    has_rigid = false
    has_disordered = false
    
    for dom in domains
        d_start, d_end = map_domain_to_dna(dom, start_codon, strand)
        # Check if the domain physically intersects the overlap zone
        if max(d_start, overlap_start) <= min(d_end, overlap_end)
            if dom.color == :midnightblue
                has_rigid = true
            elseif dom.color == :cornflowerblue
                has_disordered = true
            end
        end
    end
    
    if has_rigid return "Rigid" end
    if has_disordered return "Disordered" end
    return "Unannotated"
end

# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
survivors_df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)
base_ncbi_dir = raw"D:\ncbi_downloads\target_genomes\survivor_data\ncbi_dataset\data"

# Trackers
results = Dict{String, Int}(
    "Rigid vs Rigid" => 0,
    "Rigid vs Disordered" => 0,
    "Rigid vs Unannotated" => 0,
    "Disordered vs Disordered" => 0,
    "Disordered vs Unannotated" => 0,
    "Unannotated vs Unannotated" => 0,
    "False Overlap (Intergenic Gap)" => 0
)

println("Scanning $(nrow(survivors_df)) gene pairs...\n")

for row in eachrow(survivors_df)
    genome_acc = row.genome
    clean_gene_A = replace(row.gene_A, "gene-" => "")
    clean_gene_B = replace(row.gene_B, "gene-" => "")
    
    genome_dir = joinpath(base_ncbi_dir, genome_acc)
    gff_files = filter(f -> endswith(f, ".gff") || endswith(f, ".gff3"), readdir(genome_dir))
    if isempty(gff_files) continue end
    gff_path = joinpath(genome_dir, gff_files[1])
    
    coords_A = get_gene_coords(gff_path, clean_gene_A)
    coords_B = get_gene_coords(gff_path, clean_gene_B)
    if coords_A === nothing || coords_B === nothing continue end
    
    # Define physical overlap
    overlap_start = max(coords_A[1], coords_B[1])
    overlap_end = min(coords_A[2], coords_B[2])
    
    # Catch non-overlapping gaps (like STM3908 vs STM3909)
    if overlap_start > overlap_end
        results["False Overlap (Intergenic Gap)"] += 1
        continue
    end
    
    # Fetch Domains
    domains_A = get(domain_cache, clean_gene_A, PlotDomain[])
    domains_B = get(domain_cache, clean_gene_B, PlotDomain[])
    
    start_codon_A = coords_A[3] == '+' ? coords_A[1] : coords_A[2]
    start_codon_B = coords_B[3] == '+' ? coords_B[1] : coords_B[2]
    
    state_A = classify_overlap_zone(domains_A, start_codon_A, coords_A[3], overlap_start, overlap_end)
    state_B = classify_overlap_zone(domains_B, start_codon_B, coords_B[3], overlap_start, overlap_end)
    
    # Sort the states so "Rigid vs Disordered" is the same as "Disordered vs Rigid"
    pair_states = sort([state_A, state_B])
    category_key = "$(pair_states[2]) vs $(pair_states[1])" # Hack to match standard dictionary keys
    
    # Clean up naming
    if category_key == "Rigid vs Rigid" results["Rigid vs Rigid"] += 1
    elseif category_key == "Rigid vs Disordered" results["Rigid vs Disordered"] += 1
    elseif category_key == "Unannotated vs Rigid" results["Rigid vs Unannotated"] += 1
    elseif category_key == "Disordered vs Disordered" results["Disordered vs Disordered"] += 1
    elseif category_key == "Unannotated vs Disordered" results["Disordered vs Unannotated"] += 1
    elseif category_key == "Unannotated vs Unannotated" results["Unannotated vs Unannotated"] += 1
    end
end

# Create Summary DataFrame
summary_df = DataFrame(
    Overlap_Structural_Status = collect(keys(results)),
    Pair_Count = collect(values(results))
)
sort!(summary_df, :Pair_Count, rev=true)

CSV.write("overlap_statistics.csv", summary_df)

println("--------------------------------------------------")
println("🏆 FINAL STRUCTURAL TUG-OF-WAR STATISTICS")
println("--------------------------------------------------")
for row in eachrow(summary_df)
    println(rpad(row.Overlap_Structural_Status, 32), " | ", row.Pair_Count)
end
println("--------------------------------------------------")
println("Exported to 'overlap_statistics.csv'")



# next, let's hunt for the mythical "Rigid vs Rigid" pair that overlaps in the same physical space. This is the unicorn of gene overlaps and would be a fascinating find if it exists!


# 1. Load the Cache

domain_cache = Dict{String, Any}()
for (k, v) in raw_cache
    domain_cache[String(k)] = v
end

function map_domain_to_dna(domain, start_codon::Int, strand::Char)
    if strand == '+'
        return (start_codon + (domain.start_aa - 1) * 3, start_codon + (domain.end_aa * 3) - 1)
    else
        return (start_codon - (domain.end_aa * 3) + 1, start_codon - (domain.start_aa - 1) * 3)
    end
end

function get_gene_coords(gff_path, gene_id)
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) ∈ ("CDS", "gene")
            if occursin(gene_id, string(GFF3.attributes(record)))
                start = GFF3.seqstart(record)
                end_pos = GFF3.seqend(record)
                strand = GFF3.strand(record) == GenomicFeatures.Strand('+') ? '+' : '-'
                close(reader)
                return (start, end_pos, strand)
            end
        end
    end
    close(reader)
    return nothing
end


println("Hunting for the impossible 'Rigid vs Rigid' pair...\n")

for row in eachrow(survivors_df)
    clean_gene_A = replace(row.gene_A, "gene-" => "")
    clean_gene_B = replace(row.gene_B, "gene-" => "")
    
    domains_A = get(domain_cache, clean_gene_A, [])
    domains_B = get(domain_cache, clean_gene_B, [])
    
    # Quick filter: If both don't have at least one rigid domain, skip it
    if !any(d -> d.color == "midnightblue", domains_A) || !any(d -> d.color == "midnightblue", domains_B)
        continue
    end
    
    genome_dir = joinpath(base_ncbi_dir, row.genome)
    gff_files = filter(f -> endswith(f, ".gff") || endswith(f, ".gff3"), readdir(genome_dir))
    if isempty(gff_files) continue end
    
    coords_A = get_gene_coords(joinpath(genome_dir, gff_files[1]), clean_gene_A)
    coords_B = get_gene_coords(joinpath(genome_dir, gff_files[1]), clean_gene_B)
    if coords_A === nothing || coords_B === nothing continue end
    
    overlap_start = max(coords_A[1], coords_B[1])
    overlap_end = min(coords_A[2], coords_B[2])
    if overlap_start > overlap_end continue end
    
    start_A = coords_A[3] == '+' ? coords_A[1] : coords_A[2]
    start_B = coords_B[3] == '+' ? coords_B[1] : coords_B[2]
    
    rigid_A_names = String[]
    for dom in domains_A
        if dom.color == "midnightblue"
            d_start, d_end = map_domain_to_dna(dom, start_A, coords_A[3])
            if max(d_start, overlap_start) <= min(d_end, overlap_end)
                push!(rigid_A_names, dom.name)
            end
        end
    end
    
    rigid_B_names = String[]
    for dom in domains_B
        if dom.color == "midnightblue"
            d_start, d_end = map_domain_to_dna(dom, start_B, coords_B[3])
            if max(d_start, overlap_start) <= min(d_end, overlap_end)
                push!(rigid_B_names, dom.name)
            end
        end
    end
    
    if !isempty(rigid_A_names) && !isempty(rigid_B_names)
        println("🦄 FOUND IT!")
        println("Genome:  ", row.genome)
        println("Gene A:  ", clean_gene_A, " -> Overlapping Domains: ", join(unique(rigid_A_names), ", "))
        println("Gene B:  ", clean_gene_B, " -> Overlapping Domains: ", join(unique(rigid_B_names), ", "))
        println("Overlap: ", row.overlap_length, " bp")
        println("-"^50)
    end
end