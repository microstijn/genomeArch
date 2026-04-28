# Setup Environment 
using Pkg
# Activates the project environment in the current directory (where Project.toml is)
project_dir = joinpath(@__DIR__, "..")
Pkg.activate(project_dir)

using CSV
using DataFrames
using CairoMakie

println("Initializing algorithmic bias analysis script...")

# ==========================================================================
# 1. DATA LOADING & FILTERING
# ==========================================================================
println("Loading overlap data...")
df = CSV.File(raw"D:\pipeline_output\overlapping_pairs_annotated.csv") |> DataFrame
dfs = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)


# Filter for opposite-strand overlaps (Convergent or Divergent)
opp_strand_df = filter(row -> row.overlap_type in ["Convergent", "Divergent"], df)

# Identify genomes with >120bp overlaps and subset the data
println("Filtering for genomes with >120bp overlaps...")
target_genomes = Set(filter(:overlap_length => >(120), opp_strand_df).genome)
subset_df = filter(:genome => g -> g in target_genomes, opp_strand_df)

# ==========================================================================
# 2. ISOLATE EXEMPLARS VS. GENERAL POPULATION
# ==========================================================================
println("Isolating exemplar genomes...")
# Identify the ~10 genomes that have experimentally verified overlaps > 120bp
exemplar_genomes = Set(filter(:overlap_length => >(120), opp_strand_df).genome)
exemplar_df = filter(:genome => g -> g in exemplar_genomes, opp_strand_df)

# Map total genes
gene_counts = Dict(dfs.accession .=> dfs.p_gene_nr .+ dfs.n_gene_nr)

# Calculate pool sizes
total_genes_pop = sum(get(gene_counts, g, 0) for g in unique(opp_strand_df.genome))
total_genes_ex = sum(get(gene_counts, g, 0) for g in exemplar_genomes)

println("Found $(length(exemplar_genomes)) exemplar genomes.")
println("Total genes in Population: $total_genes_pop | Total genes in Exemplars: $total_genes_ex")

# ==========================================================================
# 3. CALCULATE EXACT FREQUENCIES (PER 1000 GENES)
# ==========================================================================
max_x = 200
lengths = collect(1:max_x)

pop_freq = zeros(Float64, max_x)
ex_freq = zeros(Float64, max_x)

for x in 1:max_x
    # Population counts
    c_pop = nrow(filter(:overlap_length => ==(x), opp_strand_df))
    pop_freq[x] = total_genes_pop > 0 ? (c_pop / total_genes_pop) * 1000.0 : 0.0
    
    # Exemplar counts
    c_ex = nrow(filter(:overlap_length => ==(x), exemplar_df))
    ex_freq[x] = total_genes_ex > 0 ? (c_ex / total_genes_ex) * 1000.0 : 0.0
end

# ==========================================================================
# 4. CALCULATE MISSING AREA (>120bp)
# ==========================================================================
missing_density_per_1000 = 0.0

for x in 120:max_x
    diff = ex_freq[x] - pop_freq[x]
    if diff > 0
        missing_density_per_1000 += diff
    end
end

total_estimated_missing = (missing_density_per_1000 / 1000.0) * total_genes_pop

println("\n*** EMPIRICAL ESTIMATION RESULTS ***")
println("Missing Density (120-200bp): $(round(missing_density_per_1000, digits=3)) per 1000 genes.")
println("Estimated Total Suppressed Overlaps across 20k genomes: ~$(round(Int, total_estimated_missing))")
println("************************************\n")

# ==========================================================================
# 5. VISUALIZATION (MONOCHROME + HIGHLIGHT)
# ==========================================================================
println("Generating plot...")
begin

set_theme!(theme_minimal(), font="Helvetica")

fig = Figure(size = (650, 450)) 
ax = Axis(fig[1, 1],
    title = "Empirical Estimation of Falsely Suppressed Genes (>120bp)",
    xlabel = "Overlap Length (bp)",
    ylabel = "Frequency per 1000 Genes",
    xgridvisible = false,
    ygridvisible = true,
)

# Background Zone
vspan!(ax, 120, max_x, color = (:grey80, 0.3), label = "Algorithmically Suppressed Zone")
vlines!(ax, [120], color = :black, linestyle = :dash, linewidth = 2)

# Shaded Area (The "Missing" Data)
band_x = 120:max_x
# Only shade where exemplar is higher than population to avoid messy visuals if noise drops it below
band_lower = pop_freq[120:max_x]
band_upper = max.(pop_freq[120:max_x], ex_freq[120:max_x])
band!(ax, band_x, band_lower, band_upper, color = (:black, 0.25), label = "Missing Overlap Area")

# Line 1: General Population (Observed)
lines!(ax, lengths, pop_freq, color = :grey50, linewidth = 2.0, label = "All Genomes (PGAP Filtered)")

# Line 2: Exemplars (Expected Baseline)
# Note: With only ~10 genomes, 1bp resolution might look jagged. 
# The band handles this, but we use a distinct black line.
lines!(ax, lengths, ex_freq, color = :black, linewidth = 2.5, label = "Exemplar Genomes (Unfiltered)")
scatter!(ax, lengths[120:end], ex_freq[120:end], color = :black, markersize = 4)

# Zoom the Y-axis to focus on the tail end. We ignore the massive 0-20bp spike for readability.
max_y_tail = maximum(ex_freq[60:max_x])
ylims!(ax, 0, max_y_tail * 10)
xlims!(ax, 0, max_x)

axislegend(ax, position = :rt)

display(fig)
end
save("empirical_missing_area_estimation.png", fig, px_per_unit = 3)
save("empirical_missing_area_estimation.pdf", fig)

println("Done! Empirical area plot saved.")