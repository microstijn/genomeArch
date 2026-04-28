using LinearAlgebra
using CairoMakie

# ==========================================
# 1. Kneedle Algorithm Implementation
# ==========================================
function find_kneedle_horizon(df_bin)
    # 1. Extract X (gap size) and Y (log overlaps)
    # Ensure data is sorted by gap size decreasing (400 down to 0) to match the curve shape
    df_sorted = sort(df_bin, :gap_bin, rev=true)
    x = Float64.(df_sorted.gap_bin)
    y = Float64.(df_sorted.log_overlap_median)
    
    # 2. Normalize X and Y to [0, 1] to prevent scale distortion
    # Because X goes from 400 to 0, we normalize so the "start" (relaxed) is 0 and "end" (compact) is 1
    x_norm = (x[1] .- x) ./ (x[1] - x[end]) 
    y_norm = (y .- minimum(y)) ./ (maximum(y) - minimum(y))
    
    # 3. Define the Secant Line (from first normalized point to last)
    p1 = [x_norm[1], y_norm[1]]
    p2 = [x_norm[end], y_norm[end]]
    
    # Precompute line length for the distance formula
    line_vec = p2 - p1
    line_length = norm(line_vec)
    
    # 4. Calculate perpendicular distance from every point to the secant line
    distances = Float64[]
    for i in 1:length(x_norm)
        p0 = [x_norm[i], y_norm[i]]
        # 2D cross product equivalent for distance from point to line
        dist = abs((p2[1] - p1[1]) * (p1[2] - p0[2]) - (p1[1] - p0[1]) * (p2[2] - p1[2])) / line_length
        push!(distances, dist)
    end
    
    # 5. Find the index of the maximum distance (The Knee)
    knee_idx = argmax(distances)
    knee_bp = x[knee_idx]
    
    return knee_bp, x_norm, y_norm, distances, knee_idx
end

println("Calculating Kneedle Horizon...")
kneedle_bp, x_norm, y_norm, dists, k_idx = find_kneedle_horizon(df_binned)

println("-------------------------------------------------")
println("📐 KNEEDLE COMPACTION HORIZON: ", kneedle_bp, " bp")
println("-------------------------------------------------")

# ==========================================
# 2. Kneedle Visualization (Supplementary Figure)
# ==========================================
fig = Figure(size = (900, 450), font = "Helvetica")

# --- Left Panel: The Normalized Geometric View ---
ax1 = Axis(fig[1, 1],
    title = "Kneedle Algorithm (Normalized Geometric View)",
    xlabel = "Normalized Genome Compaction (0 = Relaxed, 1 = Max Compaction)",
    ylabel = "Normalized Log10(Overlaps)",
    titlesize = 16
)

# Plot Secant Line
lines!(ax1, [x_norm[1], x_norm[end]], [y_norm[1], y_norm[end]], color = :gray50, linestyle = :dash, linewidth = 2, label = "Secant Line")

# Plot Normalized Data Curve
lines!(ax1, x_norm, y_norm, color = :dodgerblue, linewidth = 3, label = "Smoothed Trajectory")
scatter!(ax1, x_norm, y_norm, color = :black, markersize = 6)

# Highlight Maximum Perpendicular Distance
p1_norm = [x_norm[1], y_norm[1]]
p2_norm = [x_norm[end], y_norm[end]]
p0_norm = [x_norm[k_idx], y_norm[k_idx]]

# Calculate projection point on the secant line for the drawing
t = dot(p0_norm - p1_norm, p2_norm - p1_norm) / norm(p2_norm - p1_norm)^2
proj = p1_norm + t * (p2_norm - p1_norm)

lines!(ax1, [p0_norm[1], proj[1]], [p0_norm[2], proj[2]], color = :red, linewidth = 3, label = "Max Distance (Knee)")
scatter!(ax1, [x_norm[k_idx]], [y_norm[k_idx]], color = :gold, marker = :star5, markersize = 20, strokewidth=1, label="Kneedle Optimum")

axislegend(ax1, position = :lt, framevisible = false)


# --- Right Panel: Mapped Back to Biological Data ---
ax2 = Axis(fig[1, 2],
    title = "Kneedle Horizon on Biological Data",
    xlabel = "Family Median Absolute Gap Mean (bp) [Decreasing →]",
    ylabel = "Log10 (Convergent overlaps + 1)",
    xreversed = true,
    titlesize = 16
)

# Raw binned data
scatter!(ax2, df_binned.gap_bin, df_binned.log_overlap_median, color = :darkorange, markersize = 12, strokewidth = 1, strokecolor = :black, label = "10-bp Binned Medians")

# Kneedle Threshold
vlines!(ax2, [kneedle_bp], color = :black, linestyle = :dash, linewidth = 2, label = "Kneedle Horizon ($kneedle_bp bp)")

# Highlight the specific bin chosen
kneedle_y = filter(r -> r.gap_bin == kneedle_bp, df_binned).log_overlap_median[1]
scatter!(ax2, [kneedle_bp], [kneedle_y], color = :gold, marker = :star5, markersize = 25, strokewidth=1, strokecolor=:black)

axislegend(ax2, position = :lt, framevisible = false)

display(fig)
save("kneedle_verification.png", fig, px_per_unit = 4)