# test/runtests.jl

# --- Setup ---
using Test
using Pkg

# Activate the project environment (the parent directory of this test folder)
project_dir = joinpath(@__DIR__, "..")
Pkg.activate(project_dir)

# Import the code to be tested from your main module
using genArch
using DataFrames
using CSV

# --- Test Data ---
# The GFF content is now a constant within the test module.
const GFF_CONTENT = """
##gff-version 3
##sequence-region test_contig 1 4000
test_contig	TEST	region	1	4000	.	.	.	ID=test_contig
test_contig	TEST	gene	100	200	.	+	.	ID=p_gene1
test_contig	TEST	gene	150	250	.	+	.	ID=p_gene2
test_contig	TEST	gene	220	320	.	+	.	ID=p_gene2a
test_contig	TEST	gene	400	500	.	-	.	ID=n_gene1
test_contig	TEST	gene	480	580	.	-	.	ID=n_gene2
test_contig	TEST	gene	550	650	.	-	.	ID=n_gene2a
test_contig	TEST	gene	700	800	.	+	.	ID=p_gene3
test_contig	TEST	gene	780	880	.	-	.	ID=n_gene3
test_contig	TEST	gene	900	1000	.	+	.	ID=p_gene3a
test_contig	TEST	gene	950	1050	.	-	.	ID=n_gene3a
test_contig	TEST	gene	1000	1100	.	-	.	ID=n_gene4
test_contig	TEST	gene	1090	1190	.	+	.	ID=p_gene4
test_contig	TEST	gene	1200	1300	.	-	.	ID=n_gene4a
test_contig	TEST	gene	1280	1380	.	+	.	ID=p_gene4a
test_contig	TEST	gene	1400	1500	.	+	.	ID=p_gene5
test_contig	TEST	gene	1700	1800	.	-	.	ID=n_gene5
test_contig	TEST	gene	2000	2100	.	+	.	ID=p_gene6
test_contig	TEST	gene	2300	2400	.	-	.	ID=n_gene6
test_contig	TEST	gene	2600	2700	.	+	.	ID=p_gene7
test_contig	TEST	gene	2900	3000	.	-	.	ID=n_gene7
##sequence-region no_genes_contig 1 1000
no_genes_contig	TEST	region	1	1000	.	.	.	ID=no_genes_contig
##sequence-region no_overlaps_contig 1 5000
no_overlaps_contig	TEST	region	1	5000	.	.	.	ID=no_overlaps_contig
no_overlaps_contig	TEST	gene	1000	1100	.	+	.	ID=p_no_overlap1
no_overlaps_contig	TEST	gene	2000	2100	.	-	.	ID=n_no_overlap1
no_overlaps_contig	TEST	gene	3000	3100	.	+	.	ID=p_no_overlap2
no_overlaps_contig	TEST	gene	4000	4100	.	-	.	ID=n_no_overlap2
"""

# --- Main Test Set ---
@testset "genArch.jl Tests" begin

    # Create temporary files and directories for the test run
    temp_gff_dir = mktempdir()
    genome_subdir = joinpath(temp_gff_dir, "test_genome")
    mkpath(genome_subdir)
    temp_gff_file = joinpath(genome_subdir, "genomic.gff")
    temp_output_file = joinpath(temp_gff_dir, "results.csv")

    try
        write(temp_gff_file, GFF_CONTENT)
        
        # --- Run the function to be tested ---
        calculate_architecture(temp_gff_dir, temp_output_file)
        
        @test isfile(temp_output_file)
        results_df = CSV.File(temp_output_file) |> DataFrame
        @test nrow(results_df) == 3

        # --- Test Set 1: Complex Contig ---
        @testset "Complex Contig Validation" begin
            row = filter(r -> r.contig_name == "test_contig", results_df)
            @test nrow(row) == 1
            
            expected = Dict(
                :p_gene_nr => 10, :n_gene_nr => 10, :p_gene_length_sum => 1000, :n_gene_length_sum => 1000,
                :p_U_overlap_nr => 2, :p_U_overlap_length_sum => 80, :n_U_overlap_nr => 3, :n_U_overlap_length_sum => 100,
                :C_overlap_nr => 2, :C_length_sum => 70, :D_overlap_nr => 2, :D_length_sum => 30,
                :p_gap_length_sum => 1680, :n_gap_length_sum => 1700, :p_gap_mean => 240.0, :p_gap_median => 100.0,
                :n_gap_mean => 1700 / 6, :n_gap_median => 265.0, :operon_nr => 4, :operonicity_score => 50.0,
                :mean_operon_size => 2.5, :strand_asymmetry => 0.5, :divergent_pairs_nr => 4, :convergent_pairs_nr => 5,
                :gene_density_gradient_std => 0.0
            )
            for (metric, val) in expected
                typeof(val) <: AbstractFloat ? (@test row[1, metric] ≈ val atol=1e-2) : (@test row[1, metric] == val)
            end
        end
        
        # --- Test Set 2: Contig with No Genes ---
        @testset "No Genes Contig Validation" begin
            row = filter(r -> r.contig_name == "no_genes_contig", results_df)
            @test nrow(row) == 1
            
            for col in names(row)
                if occursin("gene", col) || occursin("overlap", col) || occursin("gap", col) || occursin("operon", col) || occursin("pairs", col)
                    @test row[1, col] == 0
                end
            end
            @test row[1, :strand_asymmetry] == 0.5
        end
        
        # --- Test Set 3: Contig with No Overlaps ---
        @testset "No Overlaps Contig Validation" begin
            row = filter(r -> r.contig_name == "no_overlaps_contig", results_df)
            @test nrow(row) == 1
            
            @test row[1, :p_gene_nr] == 2
            @test row[1, :n_gene_nr] == 2
            @test row[1, :p_gene_length_sum] == 200
            @test row[1, :n_gene_length_sum] == 200
            
            for col in names(row)
                if occursin("overlap", col) || occursin("operon", col) || occursin("pairs", col)
                    @test row[1, col] == 0
                end
            end
            @test row[1, :p_gap_length_sum] == 1900
            @test row[1, :n_gap_length_sum] == 1900
        end

    finally
        # Cleanup the temporary directory and files
        rm(temp_gff_dir, recursive=true, force=true)
    end
end
