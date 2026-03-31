#=====================================================
# Description:   Example script to run the genArch package pipeline.
# Author:        SHP
# Date:          2025
=====================================================#

# Setup Environment 
using Pkg
# Activates the project environment in the current directory (where Project.toml is)
project_dir = joinpath(@__DIR__, "..")
Pkg.activate(project_dir)

# Import genAarch
# This line gives you access to the exported functions from your modules.
using Revise
using genArch

# Define I/O
data_dir = raw"D:\ncbi_downloads\bactera_reference\ncbi_dataset\data"
taxdump_dir = "D:/ncbi_downloads/taxdump/"
output_dir = raw"D:\pipeline_output"

# Create the output directory if it doesn't exist
mkpath(output_dir)

# Process taxonomy
jsonl_files = [joinpath(data_dir, "assembly_data_report.jsonl")]
process_taxonomy(jsonl_files, taxdump_dir, output_dir)

# Fetch environments
tax_output_file = joinpath(output_dir, "assembly_data_report_TaxId.csv")
fetch_environments(tax_output_file, output_dir)

 #Calculate architecture
gff_dir = data_dir
arch_output_file = joinpath(output_dir, "genome_architecture_metrics.csv")
calculate_architecture(gff_dir, arch_output_file)

consolidate_to_genomes(
    arch_output_file,
    joinpath(output_dir, "per_genome_architecture_metrics.csv")
)

calculate_density_score(
    joinpath(output_dir, "per_genome_architecture_metrics.csv"),
    joinpath(output_dir, "density_scores.csv")
)

perform_pic_analysis(
    joinpath(output_dir, "density_scores.csv"),
    raw"D:\pipeline_output\assembly_data_report_TaxId.csv",
    raw"D:\pipeline_output\assembly_genome_environments.tsv"
)



using CSV
using DataFrames
f = CSV.File(joinpath(output_dir, "per_genome_architecture_metrics.csv")) |> DataFrame


prune_gtdb_tree(
    raw"D:\GTDB\ar53_r220.tree",
    joinpath(output_dir, "per_genome_architecture_metrics.csv"),
    
)

inspect_tree_file(
    raw"D:\GTDB\ar53_r220.tree"
)

tree = open(parsenewick, Phylo.path(raw"D:\GTDB\ar53_r220.tree"))

using PhyloNetworks
net = readnewick(readlines(raw"D:\GTDB\ar53_r220.tree"));
readstring(raw"D:\GTDB\ar53_r220.tree")

using NewickTree
tree_string = read(raw"D:\GTDB\ar53_r220.tree", String)

readTopology(tree_string)

run_all_tests()



# ------------------------- 
# test the new functionality. 

tempura = raw"C:\Users\peete074\OneDrive - Wageningen University & Research\programming\genomeArch\tempura\200617_TEMPURA.csv"
using CSV
using DataFrames
df_tempura = CSV.File(tempura, quoted = false) |> DataFrame

df_tempura = CSV.File(tempura, quoted =false, silencewarnings=true) |> DataFrame


# Lad both datasets
genarch = joinpath(output_dir, "per_genome_architecture_metrics.csv")
tax = joinpath(output_dir, "assembly_data_report_TaxId.csv")

df_genarch = CSV.File(genarch) |> DataFrame
df_tax = CSV.File(tax) |> DataFrame

rename!(df_genarch, :genome_name => :accession)

# Keep only the columns you need from the tax report to avoid bloat
df_tax_subset = select(df_tax, :accession, :taxId) 

# Merge them together (assuming 'accession' is the shared column name in both)
df_merged = leftjoin(df_genarch, df_tax_subset, on=:accession)

# Save the corrected file
CSV.write(joinpath(output_dir, "genarch_with_taxid.csv"), df_merged)

println("Successfully added taxId! Now ready for PipelineTools.jl")



df_tempura = merge_and_impute_ogt(
    joinpath(output_dir, "genarch_with_taxid.csv"),
    tempura,
    raw"D:\ncbi_downloads\taxdump\nodes.dmp",
    raw"D:\ncbi_downloads\taxdump\names.dmp",
    joinpath(output_dir, "merged_imputed_ogt.csv")
)

df_ready = merge_and_impute_lifestyle(
    df_tempura,
    raw"D:\pipeline_output\assembly_genome_environments.tsv",
    raw"D:\ncbi_downloads\taxdump\nodes.dmp",
    raw"D:\ncbi_downloads\taxdump\names.dmp"
)

best_model, processed_df = optimize_models(df_ready)
println(names(processed_df))
mechanistic_engine(processed_df)

using DataFrames, Statistics

function analyze_compression_by_lifestyle(df::DataFrame, threshold::Float64)
    println("Analyzing genomes compressed below the $(round(threshold, digits=2)) bp threshold...")
    
    # Categorize based on the imputed probability
    df.lifestyle_category = ifelse.(df.is_free_living .> 0.5, "Free-Living", "Host-Associated")
    
    # Flag genomes that have crossed the toxicity threshold
    df.is_compressed = df.mean_gap_size .< threshold
    
    # Group and summarize the statistics
    summary_df = combine(groupby(df, :lifestyle_category), 
        nrow => :Total_Genomes,
        :is_compressed => sum => :Genomes_Below_Threshold,
        :is_compressed => (x -> round(mean(x) * 100, digits=2)) => :Percent_Compressed
    )
    
    return summary_df
end

# Just pass in your DataFrame and the threshold we just found
compression_stats = analyze_compression_by_lifestyle(processed_df, 135.1182)
println(compression_stats)