using CSV
using DataFrames
using GenomicFeatures
using HTTP
using JSON3

println("Initializing JSON Domain Cache Builder...")

struct ProteinDomain
    name::String
    start_aa::Int
    end_aa::Int
end

function fetch_protein_domains(protein_id::AbstractString)
    if protein_id == "" return ProteinDomain[] end
    clean_id = replace(protein_id, r"^cds-" => "")
    url = "https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/RefSeq/$clean_id"
    domains = ProteinDomain[]
    
    try
        response = HTTP.get(url, retry=true, retries=3, readtimeout=10)
        if response.status == 200
            data = JSON3.read(String(response.body))
            if hasproperty(data, :results)
                for result in data[:results]
                    domain_name = result[:metadata][:name]
                    for prot in result[:proteins]
                        for loc in prot[:entry_protein_locations]
                            for frag in loc[:fragments]
                                push!(domains, ProteinDomain(domain_name, frag[:start], frag[:end]))
                            end
                        end
                    end
                end
            end
        end
    catch e
        # Silently catch API 404s
    end
    return domains
end



function extract_protein_ids(gff_path, gene_A_id, gene_B_id)
    prot_A, prot_B = "", ""
    reader = GFF3.Reader(open(gff_path, "r"))
    for record in reader
        if GFF3.featuretype(record) ∈ ("CDS", "gene")
            attr_str = string(GFF3.attributes(record))
            
            # UPDATED REGEX: Matches Julia's ["protein_id" => ["..."]] syntax
            m_prot = match(r"\"protein_id\"\s*=>\s*\[\"([^\"]+)\"\]", attr_str)
            extracted_prot = m_prot !== nothing ? m_prot.captures[1] : ""
            
            if occursin(gene_A_id, attr_str) && extracted_prot != ""
                prot_A = extracted_prot
            elseif occursin(gene_B_id, attr_str) && extracted_prot != ""
                prot_B = extracted_prot
            end
        end
    end
    close(reader)
    return prot_A, prot_B
end

# --- MAIN EXECUTION ---
df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)
base_ncbi_dir = raw"D:\ncbi_downloads\target_genomes\survivor_data\ncbi_dataset\data"

# Create a master dictionary to hold our JSON structure
cache_dict = Dict{String, Any}()

println("Starting API calls (this will take a few minutes)...")

for (i, row) in enumerate(eachrow(df))
    genome_acc = row.genome
    clean_gene_A = replace(row.gene_A, "gene-" => "")
    clean_gene_B = replace(row.gene_B, "gene-" => "")
    
    # Generate the Unique ID Key
    uid = "$(genome_acc)_$(clean_gene_A)_vs_$(clean_gene_B)"
    
    genome_dir = joinpath(base_ncbi_dir, genome_acc)
    if !isdir(genome_dir) continue end
    
    gff_files = filter(f -> endswith(f, ".gff"), readdir(genome_dir))
    if isempty(gff_files) continue end
    gff_path = joinpath(genome_dir, gff_files[1])
    
    prot_A, prot_B = extract_protein_ids(gff_path, clean_gene_A, clean_gene_B)
    
    print("\rProcessing $i / $(nrow(df)) : $prot_A & $prot_B ...        ")
    
    doms_A = fetch_protein_domains(prot_A)
    sleep(0.3)
    doms_B = fetch_protein_domains(prot_B)
    sleep(0.3)
    
    # Store the results under the unique ID key
    cache_dict[uid] = Dict(
        "prot_A" => prot_A,
        "prot_B" => prot_B,
        "domains_A" => [(name=d.name, start_aa=d.start_aa, end_aa=d.end_aa) for d in doms_A],
        "domains_B" => [(name=d.name, start_aa=d.start_aa, end_aa=d.end_aa) for d in doms_B]
    )
end

# Write the entire dictionary to a standalone JSON file
open("domain_cache.json", "w") do f
    JSON3.write(f, cache_dict)
end

println("\nDone! All domains independently cached to 'domain_cache.json'")