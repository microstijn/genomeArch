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
output_dir = joinpath(project_dir, "pipeline_output/")

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



using CSV
using DataFrames
f = CSV.File(raw"C:\Users\peete074\OneDrive - Wageningen University & Research\programming\genArch\pipeline_output\genome_architecture_metrics.csv") |> DataFrame

print(names(f))