using CSV
using DataFrames

df = CSV.read(raw"D:\pipeline_output\survivor_overlaps.csv", DataFrame)
# Extract all Gene A and Gene B names, strip the "gene-" prefix, and get unique values
all_genes = unique(vcat(
    replace.(df.gene_A, "gene-" => ""), 
    replace.(df.gene_B, "gene-" => "")
))

open("my_protein_list.txt", "w") do f
    for gene in all_genes
        println(f, gene)
    end
end
println("Saved $(length(all_genes)) IDs to my_protein_list.txt")