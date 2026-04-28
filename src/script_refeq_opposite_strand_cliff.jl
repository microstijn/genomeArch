
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
println("Loading data...")
df = CSV.File(raw"D:\pipeline_output\overlapping_pairs_annotated.csv") |> DataFrame
dfs = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)



# Filter for opposite-strand overlaps (Convergent or Divergent)
opp_strand_df = filter(row -> row.overlap_type in ["Convergent", "Divergent"], df)

# ==========================================================================
# 2. DATA PROCESSING (DENSITY CALCULATION)
# ==========================================================================

begin
lengths = Int[]
countss = Int[]
densities = Float64[]

println("Calculating lengths and overlaps per genome...")
for i in 1:1:200
    s_df = filter(:overlap_length => >=(i), opp_strand_df)
    pairs = nrow(s_df)
    
    # Calculate unique genomes and density
    if pairs > 0
        unique_gens = length(unique(s_df.genome))
        density = pairs / unique_gens
    else
        density = 0.0
    end
    
    push!(lengths, i)
    push!(countss, pairs)
    push!(densities, density)
end

# ==========================================================================
# 3. VISUALIZATION (TWO PANELS)
# ==========================================================================
println("Generating CairoMakie plot...")

set_theme!(theme_minimal(), font="Helvetica")

fig = Figure(size = (400, 400))

# --- PANEL 1: Total Overlapping Pairs (Top) ---
ax1 = Axis(fig[1, 1],
    title = "Algorithmic Suppression of Overlapping Genes (> 120bp)",
    ylabel = "Total Gene Pairs (Log10)",
    yscale = log10,
    xgridvisible = false,
    ygridvisible = true,
    yminorgridvisible = true,
    xticklabelsvisible = false
)

vspan!(ax1, 120, 200, color = (:crimson, 0.1), label = "Algorithmically Suppressed Zone")
vlines!(ax1, [120], color = :crimson, linestyle = :dash, linewidth = 2.5)

lines!(ax1, lengths, countss, color = :midnightblue, linewidth = 3.0)
scatter!(ax1, lengths, countss, color = :dodgerblue, markersize = 6)

ylims!(ax1, 10, 1_000_000)

# --- PANEL 2: Density / Overlaps per Genome (Bottom) ---
ax2 = Axis(fig[2, 1],
    xlabel = "Minimum Overlap Length (bp)",
    ylabel = "Mean Overlaps per Genome",
    xgridvisible = false,
    ygridvisible = true,
    yminorgridvisible = false
)

vspan!(ax2, 120, 200, color = (:crimson, 0.1))
vlines!(ax2, [120], color = :crimson, linestyle = :dash, linewidth = 2.5)

# Plot density on a linear scale
lines!(ax2, lengths, densities, color = :darkorange, linewidth = 3.0)
scatter!(ax2, lengths, densities, color = :orange, markersize = 6)

# Dynamic Y-limits based on the max density to ensure it fits cleanly
ylims!(ax2, 0, maximum(densities) * 1.2)

# --- FINAL FORMATTING ---
linkxaxes!(ax1, ax2)
xlims!(ax1, 0, 205)
rowgap!(fig.layout, 10)

display(fig)
end
save("overlap_density_multipanel.png", fig, px_per_unit = 3)
save("overlap_density_multipanel.pdf", fig)

println("Done! Multi-panel density plot saved.")