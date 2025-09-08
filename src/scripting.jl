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