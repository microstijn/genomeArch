# Filename: debug_intervaltree.jl

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
# Filename: debug_insert.jl

# Filename: debug_assignment.jl

# Filename: debug_final.jl

# 1. Import necessary packages
using GFF3
using GenomicFeatures
using IntervalTrees

# Filename: debug_working.jl
# Filename: debug_genomicfeatures.jl

# 1. Import necessary packages
using GFF3
using GenomicFeatures
# IntervalTrees is an underlying dependency of GenomicFeatures
using IntervalTrees

println("--- MWE using the high-level GenomicFeatures.jl API ---")

# Filename: debug_working.jl

# 1. Import necessary packages
using GFF3
using GenomicFeatures
# IntervalTrees is an underlying dependency of GenomicFeatures
using IntervalTrees

println("--- MWE using the high-level GenomicFeatures.jl API ---")

# 2. Fake GFF file in memory
fake_gff_content = """
##gff-version 3
##sequence-region contig1 1 5000
contig1	TEST	gene	500	1000	.	+	.	ID=gene_pos_1
contig1	TEST	gene	1200	1800	.	+	.	ID=gene_pos_2
contig1	TEST	gene	950	1500	.	-	.	ID=gene_neg_1
"""

# 3. Parse the fake GFF string
genes = GFF3.Record[]
reader = GFF3.Reader(IOBuffer(fake_gff_content))
for record in reader
    if GFF3.featuretype(record) == "gene"
        push!(genes, record)
    end
end
close(reader)

# 4. Extract coordinates for positive and negative strand genes
p_genes = [g for g in genes if GFF3.strand(g) == STRAND_POS]
p_starts = GFF3.seqstart.(p_genes)
p_ends = GFF3.seqend.(p_genes)

n_genes = [g for g in genes if GFF3.strand(g) == STRAND_NEG]
n_starts = GFF3.seqstart.(n_genes)
n_ends = GFF3.seqend.(n_genes)

println("\n[DEBUG] Found ", length(p_genes), " positive-strand and ", length(n_genes), " negative-strand genes.")

# 5. Test the high-level overlap detection pattern
try
    println("\n>>> Creating vectors of GenomicFeatures.Interval objects...")

    # Create vectors of GenomicFeatures.Intervals
    p_intervals = [GenomicFeatures.Interval("contig1", s, e, STRAND_POS) for (s, e) in zip(p_starts, p_ends)]
    n_intervals = [GenomicFeatures.Interval("contig1", s, e, STRAND_NEG) for (s, e) in zip(n_starts, n_ends)]

    println("\n[DEBUG] Successfully created vectors of GenomicFeatures intervals.")

    # Use the high-level eachoverlap function directly on the vectors.
    println("\n>>> Finding overlaps with `eachoverlap`...")
    
    overlap_count = 0
    # FINAL FIX: Access tuple elements with overlap[1] and overlap[2]
    for overlap in eachoverlap(p_intervals, n_intervals)
        overlap_count += 1
        println("  [Overlap Found]: ", overlap[1], " overlaps with ", overlap[2])
    end

    println("\n✅ SUCCESS: The `eachoverlap` operation completed.")
    println("[INFO] Total overlaps found: ", overlap_count)

catch e
    println("\n❌ ERROR: The operation failed. See details below.")
    showerror(stdout, e)
end

println("\n--- MWE Finished ---")