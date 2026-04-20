using CairoMakie
using DataFrames
using CSV
using Loess
using Statistics

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

df.total_genes = df.p_gene_nr .+ df.n_gene_nr
dropmissing!(df, [:absolute_gap_mean, :OGT, :coding_density_pct, :mean_gene_length])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.OGT) && row.total_genes > 100, df)

# Normalize topological counts per 1000 genes
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_C = (df.C_overlap_nr ./ df.total_genes) .* 1000.0

# Define the Temperature Cohorts
cohorts = [
    ("Psychro/Mesophiles (OGT < 30°C)", row -> row.OGT < 30.0, :steelblue),
    ("Moderate (30°C ≤ OGT < 45°C)", row -> 30.0 <= row.OGT < 45.0, :forestgreen),
    ("Thermophiles (OGT ≥ 45°C)", row -> row.OGT >= 45.0, :firebrick)
]

# ==============================================================================
# 2. SETUP THE MULTI-PANEL FIGURE
# ==============================================================================
fig = Figure(size = (1600, 1000), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Universal Compaction Mechanics Across Thermal Regimes", fontsize = 26, font = :bold)

panels = [
    (:coding_density_pct, "1. 1D Physical Density Ceiling", "Density (%)", 1, 1),
    (:mean_gene_length, "2. Protein Shrinkage (The U-Curve)", "Mean Length (bp)", 1, 2),
    (:norm_U, "3. Unidirectional Overlaps", "Events / 1k Genes", 2, 1),
    (:norm_C, "4. Convergent Collisions", "Events / 1k Genes", 2, 2)
]

axes = []
eval_points = collect(400.0:-1.0:0.0)
bin_edges = 0.0:15.0:400.0
bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2.0

# ==============================================================================
# 3. PLOT COHORT TRAJECTORIES
# ==============================================================================
for (sym, title_str, y_label, row, col) in panels
    ax = Axis(fig[row, col], 
              title = title_str, titlealign = :left,
              ylabel = y_label,
              xreversed = true, 
              limits = ((0, 400), (nothing, nothing)))
    push!(axes, ax)
    
    # Loop through each thermal cohort and plot its distinct trajectory
    for (label_str, filter_func, color) in cohorts
        sub_df = filter(filter_func, df)
        
        x_data = sub_df.absolute_gap_mean
        y_data = Float64.(sub_df[!, sym])
        
        valid_centers, bin_means = Float64[], Float64[]
        for i in 1:(length(bin_edges)-1)
            mask = bin_edges[i] .<= x_data .< bin_edges[i+1]
            if sum(mask) >= 5 # Dynamic Bounding specific to the cohort
                push!(valid_centers, bin_centers[i])
                push!(bin_means, mean(y_data[mask]))
            end
        end
        
        # Fit Loess strictly within this cohort's data bounds
        if length(valid_centers) > 3
            model = loess(valid_centers, bin_means, span=0.45)
            safe_eval = filter(x -> minimum(valid_centers) <= x <= maximum(valid_centers), eval_points)
            y_smooth = predict(model, safe_eval)
            
            # Draw the smoothed trajectory without scatter points for a clean overlay
            lines!(ax, safe_eval, y_smooth, linewidth = 5, color = color, label = label_str)
        end
    end
    
    # Add legend to the first panel only
    if row == 1 && col == 1
        axislegend(ax, position = :lt, framevisible = false, labelsize = 16)
    end
end

linkxaxes!(axes...)
for ax in axes[1:2] hidexdecorations!(ax, grid = false, ticks = false) end

axes[3].xlabel = "Evolutionary Timeline: Absolute Gap Mean (bp) [Decreasing →]"
axes[4].xlabel = "Evolutionary Timeline: Absolute Gap Mean (bp) [Decreasing →]"

rowgap!(fig.layout, 15)
colgap!(fig.layout, 20)

# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\universal_compaction_mechanics.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)