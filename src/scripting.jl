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


using CSV
using DataFrames

genarch = joinpath(output_dir, "per_genome_architecture_metrics.csv")
tax = joinpath(output_dir, "assembly_data_report_TaxId.csv")

df_genarch = CSV.File(genarch) |> DataFrame
df_tax = CSV.File(tax) |> DataFrame

rename!(df_genarch, :genome_name => :accession)

# Keep accession, taxId, and the full 7-rank taxonomic lineage
# (We intersect with actual names just in case a rank column is missing from the file)
cols_to_keep = intersect(names(df_tax), [
    "accession", "taxId", 
    "superkingdom", "phylum", "class", "order", "family", "genus", "species"
])
df_tax_subset = select(df_tax, cols_to_keep) 

# Merge them together 
df_merged = leftjoin(df_genarch, df_tax_subset, on=:accession)

# Save the corrected file
CSV.write(joinpath(output_dir, "genarch_with_full_taxonomy.csv"), df_merged)
println("Successfully added taxId! Now ready for PipelineTools.jl")

using CSV

df = CSV.File(joinpath(output_dir, "genarch_with_full_taxonomy.csv")) |> DataFrame  



df_tempura = merge_and_impute_ogt(
    joinpath(output_dir, "genarch_with_full_taxonomy.csv"),
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

using DataFrames

"""
    dereplicate_by_species(df::DataFrame)

Reduces the dataset to one representative per species to eliminate taxonomic bias.
Retains highly novel isolates (missing species) by treating them as unique.
"""
function dereplicate_by_species(df::DataFrame)
    println("Original genome count: ", nrow(df))
    
    # Create a grouping key: Use the species name if valid, otherwise fallback to the unique accession
    derep_key = [
        (ismissing(s) || s == "NA" || s == "") ? a : s 
        for (s, a) in zip(df.species, df.accession)
    ]
    
    # Add key to a temporary dataframe
    df_temp = copy(df)
    df_temp.derep_key = derep_key
    
    # Sort so that if we pick 'first', we get the most "average" genome size for that species, 
    # or you can just let it pick randomly. Sorting by genome_size is a safe default.
    sort!(df_temp, :genome_size)
    
    # Group by the key and take the middle/first representative
    # Taking the first after sorting gives the smallest, but let's just grab the first row for speed and randomness
    df_derep = combine(groupby(df_temp, :derep_key), first)
    
    # Clean up the temporary column
    select!(df_derep, Not(:derep_key))
    
    println("De-replicated genome count: ", nrow(df_derep))
    return df_derep
end

# Run it on your processed data
df_s = dereplicate_by_species(df_ready)

df_s.total_genes = df_s.p_gene_nr .+ df_s.n_gene_nr
#df_s.mean_gene_length = (df_s.p_gene_length_sum .+ df_s.n_gene_length_sum) ./ df_s.total_genes
engine_func, threshold, df = mechanistic_engine(df_s)


df = mechanistic_engine_all(df_s)
df.mean_gap_size

println(names(df_s))
df_s.total_U_overlaps = df_s.p_U_overlap_nr .+ df_s.n_U_overlap_nr
df_s.total_U_overlaps_density = (df_s.total_U_overlaps ./ df_s.total_genes) .* 1000.0

df_s.total_D_overlaps = (df_s.D_overlap_nr ./ df_s.total_genes) .* 1000.0

best_model, processed_df = optimize_models(df_s, threshold, :total_D_overlaps)


# Phase 3: The CairoMakie Master Figure (Quadratic & Logistic Edition)

using CairoMakie
using GLM
using DataFrames
using Statistics

"""
    plot_master_figure(df::DataFrame, model, threshold::Float64, output_file::String)

Generates a publication-grade 4-Panel Master Performance Figure using CairoMakie.
A: Hexbin Density Parity.
B: Binned Averages (Noise Filtered, Dynamic Binning).
C: Single-Gene Logistic Probability (The S-Curve).
D: Residuals vs. Gap Size (Flat Cluster confirming threshold fit).
"""
function plot_master_figure(df::DataFrame, model, threshold::Float64, output_file::String)
    println("Generating 4-Panel Master Figure...")

    # --- 1. Data Preparation ---
    # Get predictions and residuals based on the winning Model 6
    pred_overlaps = predict(model, df)
    actual_overlaps = df.C_overlap_density
    residuals = actual_overlaps .- pred_overlaps
    
    # Global Parity Line bounds
    min_val = min(minimum(pred_overlaps), minimum(actual_overlaps), 0.0)
    max_val = max(maximum(pred_overlaps), maximum(actual_overlaps), 110.0)
    
    # Calculate R-squared of this specific run
    current_r2 = round(r2(model), digits=4)
    println("Main Model R² = $current_r2")

    # --- 2. Figure Setup ---
    theme = Theme(fontsize=20, font="Helvetica")
    set_theme!(theme)
    fig = Figure(size = (1400, 1100))
    
    # ==========================================
    # PANEL A: The Hexbin Density Plot
    # ==========================================
    axA = Axis(fig[1, 1], 
        title = "A. Global Parity (Hexbin Density)",
        xlabel = "Predicted Overlaps (Model 6)", 
        ylabel = "Actual Empirical Overlaps",
        aspect = 1,
        xticks = [0, 20, 40, 60, 80, 100, 120],
        yticks = [0, 20, 40, 60, 80, 100, 120],
        #xscale = sqrt,
        #yscale = sqrt
        )
    
    hb = hexbin!(axA, pred_overlaps, actual_overlaps, cellsize=1.5, colormap=:inferno)
    lines!(axA, [min_val, max_val], [min_val, max_val], color=:white, linewidth=2, linestyle=:dash)
    Colorbar(fig[1, 1][1, 2], hb, label="Number of Genomes")

    # ==========================================
    # PANEL B: Binned Averages (Noise Filtered)
    # ==========================================
    axB = Axis(fig[1, 2], 
        title = "B. Binned Trend (Noise Filtered)",
        xlabel = "Predicted Overlaps (Binned)", 
        ylabel = "Mean Actual Overlaps",
        xticks = [10, 20, 30, 40, 50, 60])

    num_bins = 15
    min_pred = minimum(pred_overlaps)
    max_pred = maximum(pred_overlaps)
    bin_edges = range(min_pred, max_pred, length=num_bins+1)
    
    bin_centers = Float64[]
    bin_means = Float64[]
    bin_errors = Float64[]

    for i in 1:num_bins
        idx = findall(x -> bin_edges[i] <= x < bin_edges[i+1], pred_overlaps)
        if length(idx) > 3
            push!(bin_centers, (bin_edges[i] + bin_edges[i+1]) / 2)
            push!(bin_means, mean(actual_overlaps[idx]))
            push!(bin_errors, std(actual_overlaps[idx]) / sqrt(length(idx))) 
        end
    end

    lines!(axB, [min_val, max_val], [min_val, max_val], color=:gray60, linewidth=2, linestyle=:dash, label="y = x")
    errorbars!(axB, bin_centers, bin_means, bin_errors, color=:black, linewidth=2, whiskerwidth=10)
    scatter!(axB, bin_centers, bin_means, color=:royalblue4, markersize=14, strokecolor=:black, strokewidth=1)
    axislegend(axB, position=:lt)

    # ==========================================
    # PANEL C: The Logistic S-Curve (Probability)
    # ==========================================
    axC = Axis(fig[2, 1], 
        title = "C. Single-Gene Overlap Probability",
        xlabel = "Mean Intergenic Spacing (bp)", 
        ylabel = "Probability of Overlap",
        yticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    )

    # 1. Clean the temp dataframe so Makie's scatter doesn't crash
    df_temp = dropmissing(df, [:mean_gap_size, :C_overlap_nr, :total_genes])
    df_temp.overlap_prob = df_temp.C_overlap_nr ./ df_temp.total_genes
    df_temp.compression_penalty_sq = max.(0.0, threshold .- df_temp.mean_gap_size) .^ 2
    
    # Fit the binomial model dynamically for the plot
    log_mod = glm(@formula(overlap_prob ~ mean_gap_size + compression_penalty_sq), 
                  df_temp, Binomial(), LogitLink(), wts=Float64.(df_temp.total_genes))

    # 2. Create a dummy range across the gap sizes to draw the smooth S-Curve
    x_gap = collect(range(-50, 300, length=200)) # 'collect' forces it to be a standard Float array
    x_pen_sq = max.(0.0, threshold .- x_gap) .^ 2
    dummy_log = DataFrame(mean_gap_size = x_gap, compression_penalty_sq = x_pen_sq)
    
    # FIX: Force the GLM predictions to be pure Float64, stripping away the 'Missing' wrapper
    y_prob_pred = Float64.(predict(log_mod, dummy_log))

    # 3. Plot background raw probabilities and the bold predictive curve
    scatter!(
        axC,
        df_temp.mean_gap_size,
        df_temp.overlap_prob,
        color=(:teal, 0.1),
        markersize=4,
        strokewidth=0
    )
    lines!(axC, x_gap, y_prob_pred, color=:firebrick, linewidth=5, label="Logistic Fit")
    vlines!(axC, [threshold], color=:black, linewidth=2, linestyle=:dash, label="Threshold ($threshold bp)")
    
    axislegend(axC, position=:rt)
    xlims!(axC, -50, 300)
    ylims!(axC, -0.05, 0.2)

    # ==========================================
    # PANEL D: Residuals vs. Gap Size
    # ==========================================
    axD = Axis(fig[2, 2], 
        title = "D. Residuals vs. Spacing (Model 6)",
        xlabel = "Mean Intergenic Spacing (bp)", 
        ylabel = "Model Residuals (Actual - Predicted)",
        xticks = [-100, 0, 100, 200, 300, 400, 500])

    valid_idx = findall(x -> -50 <= x <= 600, df.mean_gap_size)
    
    scatter!(axD, df.mean_gap_size[valid_idx], residuals[valid_idx], color=(:midnightblue, 0.2), markersize=5, strokewidth=0)
    hlines!(axD, [0.0], color=:black, linewidth=2, linestyle=:dash)
    vlines!(axD, [threshold], color=:red2, linewidth=4, label="Threshold ($threshold bp)")
    axislegend(axD, position=:rt)

    # --- 3. Finalize and Save Output ---
    save(output_file, fig)
    println("Master Figure successfully saved to $output_file")
    
    set_theme!() # Reset theme
    return fig
end

# Run the function
plot_master_figure(processed_df, best_model, threshold, "Figure_Model6_MasterPerformance.png")


using GLM


using GLM
using DataFrames
using Statistics

println("=== Calculating the Biophysical Threshold ===")

# 1. Prepare the data: Ensure no missing values
df_logistic = dropmissing(df_s, [:mean_gap_size, :C_overlap_nr, :total_genes])
filter!(row -> !isnan(row.mean_gap_size) && row.total_genes > 0, df_logistic)

# 2. Calculate the probability (proportion) of Convergent overlaps per genome
# GLM Binomial with proportions requires the 'wts' argument (total trials/gene pairs)
df_logistic.overlap_prob = df_logistic.C_overlap_nr ./ df_logistic.total_genes

# 3. Fit the Binomial Logistic Regression (Strictly single-variable)
println("Fitting the Logistic Probability Model...")
logistic_mod = glm(
    @formula(overlap_prob ~ mean_gap_size), 
    df_logistic, 
    Binomial(), 
    LogitLink(), 
    wts=Float64.(df_logistic.total_genes)
)

# 4. View the Model Statistics
println("\n--- Model Coefficients ---")
println(coeftable(logistic_mod))

# 5. Mathematically Extract the Inflection Point (where P = 0.5)
# Equation: 0 = β0 + β1*x  =>  x = -β0 / β1
beta_0 = coef(logistic_mod)[1] # Intercept
beta_1 = coef(logistic_mod)[2] # Coefficient for mean_gap_size

threshold = -beta_0 / beta_1

println("\n--- Biophysical Threshold Result ---")
println("Intercept (β0): ", round(beta_0, digits=6))
println("Gap Size Slope (β1): ", round(beta_1, digits=6))
println("Calculated Inflection Point (P=0.5): ", round(threshold, digits=2), " bp")


# 1. Calculate the predictions and add them to your dataframe
processed_df.predicted_overlaps = predict(best_model, processed_df)

# 2. Draw a box around the "Plume" 
# (Where the model predicts < 40, but the actual reality is > 60)
df_plume = filter(row -> row.D_overlap_nr > 30, df_s)

# 3. WHO ARE THEY? (Let's check the phylogeny)
println("--- Phyla in the Underestimation Plume ---")
println(combine(groupby(df_plume, :phylum), nrow => :count))

# 4. HOW DO THEY LIVE? (Let's check their lifestyle)
println("\n--- Lifestyles in the Underestimation Plume ---")
println(combine(groupby(df_plume, :lifestyle_tier), nrow => :count))


using Statistics

println("=== Testing if OGT drives physical compression ===")
df_ogt = dropmissing(df_s, [:OGT, :mean_gap_size, :genome_size])

cor_gap = cor(df_ogt.OGT, df_ogt.mean_gap_size)
cor_size = cor(df_ogt.OGT, df_ogt.genome_size)

println("Correlation (OGT vs Mean Gap Size): ", round(cor_gap, digits=4))
println("Correlation (OGT vs Genome Size):   ", round(cor_size, digits=4))




using GLM
using CairoMakie
using Statistics

#println("=== Testing Overlap Length vs Spatial Compression ===")

# Drop missing values for the length columns and space
df_len = dropmissing(df_s, [:mean_gap_size, :C_length_sum, :n_U_overlap_length_sum, :p_U_overlap_length_sum, :D_length_sum])


df_len.C_overlap_len = df_len.C_length_sum ./ df_len.C_overlap_nr
df_len.U_overlap_len = (df_len.n_U_overlap_length_sum .+ df_len.p_U_overlap_length_sum) ./ df_len.total_U_overlaps
df_len.D_overlap_len = df_len.D_length_sum ./ df_len.D_overlap_nr

# 4. Filter out the NaNs and Infs (the division by zero artifacts)
filter!(row -> 
    !isnan(row.C_overlap_len) && !isinf(row.C_overlap_len) &&
    !isnan(row.U_overlap_len) && !isinf(row.U_overlap_len) &&
    !isnan(row.D_overlap_len) && !isinf(row.D_overlap_len) &&
    !isnan(row.mean_gap_size) && !isinf(row.mean_gap_size), 
    df_len
)

# Calculate correlations between Gap Size and Overlap Length
cor_C_len = cor(df_len.mean_gap_size, df_len.C_overlap_len)
cor_U_len = cor(df_len.mean_gap_size, df_len.U_overlap_len)
cor_D_len = cor(df_len.mean_gap_size, df_len.D_overlap_len)

println("Correlation (Space vs C Length): ", round(cor_C_len, digits=4))
println("Correlation (Space vs U Length): ", round(cor_U_len, digits=4))
println("Correlation (Space vs D Length): ", round(cor_D_len, digits=4))

# Plot the Lengths
fig = Figure(size=(1200, 400))
pub_theme = Theme(fontsize=16, font="Helvetica", Axis=(xgridvisible=false, ygridvisible=false))
set_theme!(pub_theme)

# Plot U Lengths
ax1 = Axis(fig[1, 1], title="Unidirectional Lengths", xlabel="Mean Gap Size (bp)", ylabel="Mean Overlap Length (bp)")
scatter!(ax1, df_len.mean_gap_size, df_len.U_overlap_len, color=(:dodgerblue, 0.3), markersize=5)
vlines!(ax1, [136.6], color=:red, linestyle=:dash)

# Plot C Lengths
ax2 = Axis(fig[1, 2], title="Convergent Lengths", xlabel="Mean Gap Size (bp)")
scatter!(ax2, df_len.mean_gap_size, df_len.C_overlap_len, color=(:forestgreen, 0.3), markersize=5)
vlines!(ax2, [136.6], color=:red, linestyle=:dash)

# Plot D Lengths
ax3 = Axis(fig[1, 3], title="Divergent Lengths", xlabel="Mean Gap Size (bp)")
scatter!(ax3, df_len.mean_gap_size, df_len.D_overlap_len, color=(:darkorange, 0.3), markersize=5)
vlines!(ax3, [136.6], color=:red, linestyle=:dash)

for ax in [ax1, ax2, ax3]
    ax.xscale = log10
    #ax.yscale = log10
end


fig
save("overlap_lengths.png", fig)
println("Saved plot to 'overlap_lengths.png'")

println(names(df_s))

df_s.known_lifestyle_tier