# genArch.jl: A Toolkit for genome architecture and environment
genArch.jl is a high-performance Julia package for analyzing prokaryotic genomes. The workflow extracts taxonomic and assembly metadata from NCBI datasets, retrieves associated environmental data, and calculates detailed genome architecture metrics from GFF3 files.

The tool is optimized for large-scale genomics, leveraging Julia's multi-threading capabilities to efficiently process thousands of genomes.

## Features
The toolkit provides a suite of functions to quantify the organizational principles of a genome.

## Architectural Metrics Calculated
Gene Counts & Lengths: Total number and cumulative length of genes on both the positive (+) and negative (-) strands.

### Overlap Analysis:

* Unidirectional (U_overlap): Number and total length of overlaps between adjacent genes on the same strand.

* Convergent (C_overlap): Number and total length of overlaps between genes on opposite strands oriented towards each other (--> <--).

* Divergent (D_overlap): Number and total length of overlaps between genes on opposite strands oriented away from each other (<-- -->).

* Intergenic Spacing: Total, mean, median, and standard deviation of non-coding gap lengths.

### Operon Prediction:

* operon_nr: The total number of predicted operons.

* operonicity_score: The percentage of genes located within operons.

* mean_operon_size: The average number of genes per operon.

### Local Gene Arrangements:

* divergent_pairs_nr: Counts adjacent gene pairs oriented <-- -->, often sharing a bidirectional promoter.

* convergent_pairs_nr: Counts adjacent gene pairs oriented --> <--.

### Chromosome-Scale Organization:

* strand_asymmetry: The ratio of genes on the positive strand to the total number of genes.

* gene_density_gradient_std: Measures the variability of gene density across the genome.

## Data Sources
All data processed by this toolkit is sourced from the National Center for Biotechnology Information (NCBI).

### Genome Assemblies and Metadata
The required GFF3 annotation files and assembly metadata can be downloaded using the NCBI datasets command-line tool.

For example, to download all available reference genomes for Bacteria (taxon: 2):

```bash
# Download the dehydrated package for Bacteria
datasets download genome taxon 2 --assembly-source refseq --reference --include gff3,assembly-report --dehydrated --filename bacteria_reference.zip

# Rehydrate the package to retrieve the data files
datasets rehydrate --directory bacteria_reference/
```

### Taxonomy Database
Full taxonomic lineage information is derived from the NCBI Taxonomy database dump files (nodes.dmp and names.dmp). These files can be obtained from the NCBI FTP server:

```bash
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```

## Installation
Clone the repository:

```bash
git clone https://github.com/microstijn/genomeArch/
cd genArch
```

### Install Dependencies:
From within the Julia REPL, activate the project environment and instantiate the required packages.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Pipeline Overview
The toolkit is designed to be run as a sequential pipeline using three main functions exported by the genArch module:

* process_taxonomy(): Parses NCBI assembly_data_report.jsonl files and a taxdump directory to generate a master file linking each genome assembly to its full taxonomic lineage.

* fetch_environments(): Takes the output from the taxonomy step and enriches it with environmental data by querying external databases for each taxon.

* calculate_architecture(): Processes a directory of GFF3 files to calculate all genome architecture metrics.

## Usage
The analysis can be run from a control script that calls the main pipeline functions in order.

### Directory Structure
The tool expects GFF files to be organized into subdirectories, where each subdirectory represents a single genome.

```
/your_gff_directory/
├── /genome_A/
│   └── genomic.gff
└── /genome_B/
    └── genomic.gff
```

### Example Workflow Script
Modify a script like the one below to point to your data directories, then execute it.

```julia
# Setup Environment 
using Pkg
project_dir = @__DIR__ # Assumes script is in the project's root
Pkg.activate(project_dir)

# Import the package modules
using Revise
using genArch

# Define I/O paths
data_dir = "/path/to/your/ncbi_dataset/data"
taxdump_dir = "/path/to/your/taxdump/"
output_dir = joinpath(project_dir, "pipeline_output/")
mkpath(output_dir)

# 1. Process taxonomy
jsonl_files = [joinpath(data_dir, "assembly_data_report.jsonl")]
process_taxonomy(jsonl_files, taxdump_dir, output_dir)

# 2. Fetch environments
tax_output_file = joinpath(output_dir, "assembly_data_report_TaxId.csv")
fetch_environments(tax_output_file, output_dir)

# 3. Calculate architecture
gff_dir = data_dir
arch_output_file = joinpath(output_dir, "genome_architecture_metrics.csv")
calculate_architecture(gff_dir, arch_output_file)
```
