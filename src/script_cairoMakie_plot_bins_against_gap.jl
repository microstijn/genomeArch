using CairoMakie
using DataFrames
using Statistics

println("Generating Final 4-Panel Manuscript Figure: The Single Engine of Compaction...")

# ---------------------------------------------------------
# 1. PREPARE THE DATA 
# ---------------------------------------------------------
df.mean_gene_length = (df.p_gene_length_sum .+ df.n_gene_length_sum) ./ df.total_genes
df_plot = dropmissing(df, [:mean_gap_size, :total_genes, :genome_size, :mean_gene_length,
                           :C_overlap_nr, :C_length_sum, 
                           :D_overlap_nr, :D_length_sum, 
                           :p_U_overlap_nr, :n_U_overlap_nr,
                           :p_U_overlap_length_sum, :n_U_overlap_length_sum])
filter!(r -> !isnan(r.mean_gap_size) && r.total_genes > 0 && r.genome_size > 0, df_plot)

# --- Metrics ---
df_plot.U_overlap_nr = df_plot.p_U_overlap_nr .+ df_plot.n_U_overlap_nr
df_plot.U_density = (df_plot.U_overlap_nr ./ df_plot.total_genes) .* 1000.0
df_plot.C_density = (df_plot.C_overlap_nr ./ df_plot.total_genes) .* 1000.0
df_plot.D_density = (df_plot.D_overlap_nr ./ df_plot.total_genes) .* 1000.0

df_plot.U_fraction = (df_plot.p_U_overlap_length_sum .+ df_plot.n_U_overlap_length_sum) ./ df_plot.genome_size
df_plot.C_fraction = df_plot.C_length_sum ./ df_plot.genome_size
df_plot.D_fraction = df_plot.D_length_sum ./ df_plot.genome_size

# ---------------------------------------------------------
# 2. 8-BIN LOGIC 
# ---------------------------------------------------------
function assign_bin(gap)
    if gap >= 200 return 1
    elseif gap >= 167 return 2
    elseif gap >= 134 return 3
    elseif gap >= 101 return 4
    elseif gap >= 68  return 5
    elseif gap >= 35  return 6
    elseif gap >= 15  return 7
    else return 8
    end
end

bin_labels = ["> 200", "167-200", "134-167", "101-134", "68-101", "35-68", "15-35", "< 15"]
df_plot.bin = assign_bin.(df_plot.mean_gap_size)

# Helper for melting Topology-specific data
function melt_topologies(df, cols)
    vcat([DataFrame(bin = df.bin, val = df[:, cols[i]], type = i) for i in 1:length(cols)]...)
end

df_box_density = melt_topologies(df_plot, [:U_density, :C_density, :D_density])
df_box_occupancy = melt_topologies(df_plot, [:U_fraction, :C_fraction, :D_fraction])

# Lengths (filtered for counts > 0)
u_mask, c_mask, d_mask = df_plot.U_overlap_nr .> 0, df_plot.C_overlap_nr .> 0, df_plot.D_overlap_nr .> 0
df_box_length = vcat(
    DataFrame(bin = df_plot.bin[u_mask], val = (df_plot.p_U_overlap_length_sum[u_mask] .+ df_plot.n_U_overlap_length_sum[u_mask]) ./ df_plot.U_overlap_nr[u_mask], type = 1),
    DataFrame(bin = df_plot.bin[c_mask], val = df_plot.C_length_sum[c_mask] ./ df_plot.C_overlap_nr[c_mask], type = 2),
    DataFrame(bin = df_plot.bin[d_mask], val = df_plot.D_length_sum[d_mask] ./ df_plot.D_overlap_nr[d_mask], type = 3)
)

# ---------------------------------------------------------
# 3. CAIROMAKIE 2x2 GRID
# ---------------------------------------------------------
set_theme!(theme_minimal())
update_theme!(font = "Arial", fontsize = 18, Axis = (titlesize=22, xlabelsize=18, ylabelsize=18))

fig = Figure(resolution = (1800, 1400))
colors = [:cornflowerblue, :firebrick, :forestgreen]

# PANEL A: Density
ax1 = Axis(fig[1, 1], title = "A. Overlap Density", ylabel = "Density (per 1000 genes)", xticks = (1:8, bin_labels))
boxplot!(ax1, df_box_density.bin, df_box_density.val, dodge = df_box_density.type, color = [colors[t] for t in df_box_density.type], show_outliers=false)
ylims!(ax1, 0, min(500, quantile(df_plot.U_density, 0.98)))

# PANEL B: Mean Overlap Length
ax2 = Axis(fig[1, 2], title = "B. Mean Overlap Length", ylabel = "Length (bp)", xticks = (1:8, bin_labels))
boxplot!(ax2, df_box_length.bin, df_box_length.val, dodge = df_box_length.type, color = [colors[t] for t in df_box_length.type], show_outliers=false)
ylims!(ax2, 0, 60)

# PANEL C: Genomic Occupancy
ax3 = Axis(fig[2, 1], title = "C. Total Genomic Occupancy", ylabel = "Fraction of Genome", xticks = (1:8, bin_labels))
boxplot!(ax3, df_box_occupancy.bin, df_box_occupancy.val, dodge = df_box_occupancy.type, color = [colors[t] for t in df_box_occupancy.type], show_outliers=false)
ylims!(ax3, 0, quantile(df_box_occupancy.val, 0.99))

# PANEL D: Mean Gene Length (The Negative Control)
ax4 = Axis(fig[2, 2], title = "D. Mean Coding Gene Length", ylabel = "Length (bp)", xticks = (1:8, bin_labels))
boxplot!(ax4, df_plot.bin, df_plot.mean_gene_length, color = (:grey40, 0.5), show_outliers=false, width=0.6)
# Set Y-limits to show the stability of gene length (usually around 900-1000bp)
ylims!(ax4, 700, 1100)

# Formatting adjustments
for ax in [ax1, ax2, ax3, ax4]
    vlines!(ax, 3.5, color = :black, linestyle = :dash, linewidth = 2.5)
    ax.xlabel = "Gap Size Bins (bp) [Relaxed → Squeezed]"
end

# Unified Legend
elements = [PolyElement(polycolor = c) for c in colors]
push!(elements, PolyElement(polycolor = :grey40))
Legend(fig[3, :], elements, ["Unidirectional", "Convergent", "Divergent", "Coding Genes"], 
       "Architectural Component", orientation = :horizontal, framevisible = false)

display(fig)