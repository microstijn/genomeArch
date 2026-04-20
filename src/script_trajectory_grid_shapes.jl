using CairoMakie
using DataFrames
using CSV
using Statistics
using Loess

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\genarch_with_full_taxonomy.csv", DataFrame)

dropmissing!(df, [:absolute_gap_mean, :coding_density_pct, :mean_gene_length])
filter!(row -> !isnan(row.absolute_gap_mean), df)
df.total_genes = df.p_gene_nr .+ df.n_gene_nr
filter!(row -> row.total_genes > 100, df)

# Normalize topological counts
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_C = (df.C_overlap_nr ./ df.total_genes) .* 1000.0
df.norm_D = (df.D_overlap_nr ./ df.total_genes) .* 1000.0
df.norm_abutting = (df.abutting_genes_nr ./ df.total_genes) .* 1000.0

# ==============================================================================
# 2. SETUP THE 3x2 MULTI-PANEL FIGURE
# ==============================================================================
fig = Figure(size = (1400, 1200), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Architectural Trajectories in Response to Compaction", fontsize = 28, font = :bold)

# Define the 6 panels: (Symbol, Title, Y-Axis Label, Grid Row, Grid Col, Color Index)
metrics = [
    (:strand_switch_rate, "1. Strand Legalization", "Strand Switches / Gene", 1, 1, 1),
    (:mean_gene_length, "2. Protein Shrinkage", "Mean Length (bp)", 1, 2, 2),
    (:norm_abutting, "3. Abutting Genes (0 bp)", "Events / 1k Genes", 2, 1, 3),
    (:norm_U, "4. Unidirectional Overlaps", "Events / 1k Genes", 2, 2, 4),
    (:norm_C, "5. Convergent Overlaps (Tail-to-Tail)", "Events / 1k Genes", 3, 1, 5),
    (:norm_D, "6. Divergent Overlaps (Head-to-Head)", "Events / 1k Genes", 3, 2, 6)
]

# We will evaluate lines strictly from 400 down to 0
eval_points = collect(400.0:-1.0:0.0)
bin_edges = 0.0:15.0:400.0
bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2.0

colors = Makie.wong_colors()
axes = []

# ==============================================================================
# 3. PLOT EACH TRAJECTORY
# ==============================================================================
for (sym, title_str, y_label, row, col, color_idx) in metrics
    ax = Axis(fig[row, col], 
              title = title_str, titlealign = :left,
              ylabel = y_label,
              xreversed = true, # Reverse X so time flows Left to Right
              limits = ((0, 400), (nothing, nothing))
    )
    push!(axes, ax)
    
    valid_mask = .!ismissing.(df[!, sym]) .&& .!isnan.(df[!, sym])
    x_data = df.absolute_gap_mean[valid_mask]
    y_data = Float64.(df[valid_mask, sym])
    
    # Calculate binned means to handle noise and extreme densities
    valid_centers, bin_means = Float64[], Float64[]
    for i in 1:(length(bin_edges)-1)
        mask = bin_edges[i] .<= x_data .< bin_edges[i+1]
        if sum(mask) >= 5
            push!(valid_centers, bin_centers[i])
            push!(bin_means, mean(y_data[mask]))
        end
    end
    
    # Fit the mathematical smoothing line
    model = loess(valid_centers, bin_means, span=0.45)
    
    # Restrict predictions to the data's bounding box to prevent Loess extrapolation errors
    safe_eval = filter(x -> minimum(valid_centers) <= x <= maximum(valid_centers), eval_points)
    y_smooth = predict(model, safe_eval)
    
    # Plot the raw binned data points (faded)
    scatter!(ax, valid_centers, bin_means, color = (:gray60, 0.6), markersize = 12)
    
    # Overlay the smooth trajectory curve (thick line)
    lines!(ax, safe_eval, y_smooth, linewidth = 5, color = colors[color_idx])
end

# Link all X-axes so they perfectly align vertically
linkxaxes!(axes...)

# Hide the X-axis labels for the top two rows to create a clean, unified dashboard
for ax in axes[1:4]
    hidexdecorations!(ax, grid = false, ticks = false)
end

# Only the bottom two panels get the X-axis label
axes[5].xlabel = "Evolutionary Timeline: Mean Absolute Gap (bp) [Decreasing →]"
axes[6].xlabel = "Evolutionary Timeline: Mean Absolute Gap (bp) [Decreasing →]"

rowgap!(fig.layout, 15)
colgap!(fig.layout, 20)

# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\trajectory_grid_shapes.png"
save(output_path, fig, px_per_unit = 3)
println("Trajectory Grid saved to: $output_path")
display(fig)