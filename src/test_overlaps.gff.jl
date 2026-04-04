# =====================================================
# Description:   Validation script for genArch.ArchitectureTools
# Author:        Gemini
# Date:          2025
# =====================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using genArch
using CSV
using DataFrames
using Test

# The GFF content includes a complex contig, an empty contig, and a simple contig with no overlaps.
# The GFF content explicitly defined with \t to prevent space-formatting errors
gff_content = join([
    "##gff-version 3",
    "##sequence-region test_contig 1 4000",
    "test_contig\tTEST\tregion\t1\t4000\t.\t.\t.\tID=test_contig",
    "test_contig\tTEST\tgene\t100\t200\t.\t+\t.\tID=p_gene1",
    "test_contig\tTEST\tgene\t150\t250\t.\t+\t.\tID=p_gene2",
    "test_contig\tTEST\tgene\t220\t320\t.\t+\t.\tID=p_gene2a",
    "test_contig\tTEST\tgene\t400\t500\t.\t-\t.\tID=n_gene1",
    "test_contig\tTEST\tgene\t480\t580\t.\t-\t.\tID=n_gene2",
    "test_contig\tTEST\tgene\t550\t650\t.\t-\t.\tID=n_gene2a",
    "test_contig\tTEST\tgene\t700\t800\t.\t+\t.\tID=p_gene3",
    "test_contig\tTEST\tgene\t780\t880\t.\t-\t.\tID=n_gene3",
    "test_contig\tTEST\tgene\t900\t1000\t.\t+\t.\tID=p_gene3a",
    "test_contig\tTEST\tgene\t950\t1050\t.\t-\t.\tID=n_gene3a",
    "test_contig\tTEST\tgene\t1000\t1100\t.\t-\t.\tID=n_gene4",
    "test_contig\tTEST\tgene\t1090\t1190\t.\t+\t.\tID=p_gene4",
    "test_contig\tTEST\tgene\t1200\t1300\t.\t-\t.\tID=n_gene4a",
    "test_contig\tTEST\tgene\t1280\t1380\t.\t+\t.\tID=p_gene4a",
    "test_contig\tTEST\tgene\t1400\t1500\t.\t+\t.\tID=p_gene5",
    "test_contig\tTEST\tgene\t1700\t1800\t.\t-\t.\tID=n_gene5",
    "test_contig\tTEST\tgene\t2000\t2100\t.\t+\t.\tID=p_gene6",
    "test_contig\tTEST\tgene\t2300\t2400\t.\t-\t.\tID=n_gene6",
    "test_contig\tTEST\tgene\t2600\t2700\t.\t+\t.\tID=p_gene7",
    "test_contig\tTEST\tgene\t2900\t3000\t.\t-\t.\tID=n_gene7",
    "##sequence-region no_genes_contig 1 1000",
    "no_genes_contig\tTEST\tregion\t1\t1000\t.\t.\t.\tID=no_genes_contig",
    "##sequence-region no_overlaps_contig 1 5000",
    "no_overlaps_contig\tTEST\tregion\t1\t5000\t.\t.\t.\tID=no_overlaps_contig",
    "no_overlaps_contig\tTEST\tgene\t1000\t1100\t.\t+\t.\tID=p_no_overlap1",
    "no_overlaps_contig\tTEST\tgene\t2000\t2100\t.\t-\t.\tID=n_no_overlap1",
    "no_overlaps_contig\tTEST\tgene\t3000\t3100\t.\t+\t.\tID=p_no_overlap2",
    "no_overlaps_contig\tTEST\tgene\t4000\t4100\t.\t-\t.\tID=n_no_overlap2"
], "\n") * "\n"

# --- Test Execution ---
function run_validation()
    temp_gff_dir = mktempdir()
    genome_subdir = joinpath(temp_gff_dir, "test_genome")
    mkpath(genome_subdir)
    temp_gff_file = joinpath(genome_subdir, "genomic.gff")
    temp_output_file = joinpath(temp_gff_dir, "results.csv")

    try
        println("1. Writing temporary test GFF file...")
        write(temp_gff_file, gff_content)

        println("2. Running calculate_architecture function...")
        calculate_architecture(temp_gff_dir, temp_output_file)

        println("3. Reading and validating results...")
        if !isfile(temp_output_file)
            @error "Test failed: Output file was not created."
            return
        end

        results_df = CSV.File(temp_output_file) |> DataFrame
        
        @test nrow(results_df) == 3 # We expect one row per contig

        # --- Test Set 1: Complex Contig ---
        Test.@testset "Complex Contig Validation" begin
            row = filter(r -> r.contig_name == "test_contig", results_df)
            @test nrow(row) == 1
            
            # Base logic verified. 1-based coordinate updates applied. (100 to 200 = 101 bp)
            expected = Dict(
                :p_gene_nr => 10, :n_gene_nr => 10, 
                :p_gene_length_sum => 1010, :n_gene_length_sum => 1010, # 10 genes * 101 bp
                :strand_asymmetry => 0.5, 
                :mean_gene_length => 101.0, :std_gene_length => 0.0,
                :nested_genes_nr => 0, 
                :strand_switch_rate => 0.65 # 13 switches / 20 genes
            )
            for (metric, val) in expected
                typeof(val) <: AbstractFloat ? (@test row[1, metric] ≈ val atol=1e-2) : (@test row[1, metric] == val)
            end

            # Ensure the new biophysical columns exist and are calculated without error
            new_metrics = [:U_coupled_nr, :U_deep_nr, :C_deep_nr, :D_deep_nr, 
                           :U_in_frame_nr, :U_out_of_frame_nr, :abutting_genes_nr,
                           :absolute_gap_mean, :coding_density_pct]
            for metric in new_metrics
                @test metric in propertynames(row)
                @test !ismissing(row[1, metric])
            end
        end
        
        # --- Test Set 2: Contig with No Genes ---
        Test.@testset "No Genes Contig Validation" begin
            row = filter(r -> r.contig_name == "no_genes_contig", results_df)
            @test nrow(row) == 1
            
            # All structurally active metrics should drop to 0
            for col in names(row)
                if occursin("gene", col) && col != "genome_name" && col != "contig_name"
                    if col == "strand_asymmetry"
                        @test row[1, col] == 0.5
                    else
                        @test row[1, col] == 0
                    end
                end
            end
        end
        
        # --- Test Set 3: Contig with No Overlaps ---
        Test.@testset "No Overlaps Contig Validation" begin
            row = filter(r -> r.contig_name == "no_overlaps_contig", results_df)
            @test nrow(row) == 1
            
            @test row[1, :p_gene_nr] == 2
            @test row[1, :n_gene_nr] == 2
            @test row[1, :p_gene_length_sum] == 202 # 2 * 101
            @test row[1, :n_gene_length_sum] == 202
            
            # Ensure overlap fields perfectly zero out
            for col in names(row)
                if occursin("overlap", col) || occursin("operon", col) || occursin("pairs", col) || occursin("deep", col) || occursin("coupled", col)
                    @test row[1, col] == 0
                end
            end
            
            # 1-based coordinate gaps: 3000 - 1100 - 1 = 1899
            @test row[1, :p_gap_length_sum] == 1899
            @test row[1, :n_gap_length_sum] == 1899
        end

    catch e
        @error "An error occurred during the test run." exception=(e, catch_backtrace())
    finally
        println("4. Cleaning up temporary files...")
        rm(temp_gff_dir, recursive=true, force=true)
    end
end

# --- Run the validation ---
run_validation()
println("\nValidation finished.")