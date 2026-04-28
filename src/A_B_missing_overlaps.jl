# Setup Environment 
using Pkg

using CSV
using DataFrames
using CairoMakie

println("Initializing cumulative log-scale estimation script...")

# ==========================================================================
# 1. DATA LOADING & FILTERING
# ==========================================================================
println("Loading data...")
df = CSV.File(raw"D:\pipeline_output\overlapping_pairs_annotated.csv") |> DataFrame
dfs = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

opp_strand_df = filter(row -> row.overlap_type in ["Convergent", "Divergent"], df)

# ==========================================================================
# 2. ISOLATE EXEMPLARS VS. GENERAL POPULATION
# ==========================================================================
println("Processing baseline distributions...")
exemplar_genomes = Set(filter(:overlap_length => >(120), opp_strand_df).genome)
exemplar_df = filter(:genome => g -> g in exemplar_genomes, opp_strand_df)

gene_counts = Dict(dfs.accession .=> dfs.p_gene_nr .+ dfs.n_gene_nr)
total_genes_pop = sum(get(gene_counts, g, 0) for g in unique(opp_strand_df.genome))
total_genes_ex = sum(get(gene_counts, g, 0) for g in exemplar_genomes)

# ==========================================================================
# 3. FAST CUMULATIVE CALCULATION & PERCENTAGE
# ==========================================================================
max_x = 200
lengths = collect(1:max_x)

pop_counts = zeros(Int, max_x)
ex_counts = zeros(Int, max_x)

for x in 1:max_x
    pop_counts[x] = nrow(filter(:overlap_length => ==(x), opp_strand_df))
    ex_counts[x] = nrow(filter(:overlap_length => ==(x), exemplar_df))
end

# Reverse cumulative sum
cum_pop_counts = reverse(cumsum(reverse(pop_counts)))
cum_ex_counts = reverse(cumsum(reverse(ex_counts)))

# Normalize
cum_pop_density = (cum_pop_counts ./ total_genes_pop) .* 1000.0
cum_ex_density = (cum_ex_counts ./ total_genes_ex) .* 1000.0

# Log transform
log_pop = log10.(cum_pop_density .+ 1.0)
log_ex  = log10.(cum_ex_density .+ 1.0)

# Percentage missing
pct_missing = zeros(Float64, max_x)
for x in 1:max_x
    if cum_ex_density[x] > 0
        missing_val = ((cum_ex_density[x] - cum_pop_density[x]) / cum_ex_density[x]) * 100.0
        pct_missing[x] = max(0.0, missing_val) 
    else
        pct_missing[x] = 0.0
    end
end

# ==========================================================================
# 4. VISUALIZATION (PAPER-INSPIRED STYLE)
# ==========================================================================
println("Generating plot...")

# Define the minimalist, paper-inspired theme
paper_theme = Theme(
    font = "Helvetica",
    fontsize = 14,
    Axis = (
        xgridvisible = false,
        ygridvisible = false,          # Removed to match the paper's clean look
        topspinevisible = false,       # Open top
        rightspinevisible = false,     # Open right
        spinewidth = 1.5,              # Thicker, bolder axes
        bottomspinecolor = :black,
        leftspinecolor = :black,
        xtickwidth = 1.5,
        ytickwidth = 1.5,
    )
)
set_theme!(paper_theme)

begin
fig = Figure(size = (400, 650), backgroundcolor = :white)

# --- PANEL A: Log-Cumulative Distribution ---
ax1 = Axis(fig[1, 1],
    ylabel = "Cumulative overlaps\n(log₁₀ per 1000 genes)",
    xticklabelsvisible = false # Hide X labels since we share the axis
)

# Panel label 'a' placed at the top left, outside the axis plotting area
Label(fig[1, 1, TopLeft()], "a", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

# The dashed threshold line
vlines!(ax1, [120], color = :black, linestyle = :dash, linewidth = 1.5)

# Shaded Area
band_x = 1:max_x
band_lower = log_pop
band_upper = max.(log_pop, log_ex)
band!(ax1, band_x, band_lower, band_upper, color = (:grey70, 0.5), label = "missing overlap space")

# Lines
lines!(ax1, lengths, log_pop, color = :grey60, linewidth = 2.5, label = "all genomes (filtered)")
lines!(ax1, lengths, log_ex, color = :mediumseagreen, linewidth = 3.0, label = "exemplars (unfiltered)") # Using the paper's green for the primary line

ylims!(ax1, 0, maximum(log_ex) * 1.1)

# Clean legend with no border box, matching the paper style
axislegend(ax1, position = :rt, framevisible = false, labelsize=12)

# --- PANEL B: Percentage Missing Gradient ---
ax2 = Axis(fig[2, 1],
    xlabel = "Minimum overlap length (bp)",
    ylabel = "Estimated missing (%)"
)

# Panel label 'b'
Label(fig[2, 1, TopLeft()], "b", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

vlines!(ax2, [120], color = :black, linestyle = :dash, linewidth = 1.5)

# Plot the percentage missing
band!(ax2, lengths, zeros(max_x), pct_missing, color = (:mediumpurple, 0.3)) # Using the paper's purple for the secondary space
lines!(ax2, lengths, pct_missing, color = :mediumpurple, linewidth = 3.0)

ylims!(ax2, 0, 105) 

# --- LINK AND FORMAT ---
linkxaxes!(ax1, ax2)
xlims!(ax1, 0, max_x)
rowgap!(fig.layout, 15)

display(fig)
save("publication_styled_overlap_analysis.png", fig, px_per_unit = 4) # Higher resolution export
save("publication_styled_overlap_analysis.pdf", fig)
end
println("Done! Publication-ready figure saved.")
