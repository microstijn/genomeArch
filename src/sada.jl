using CSV
using DataFrames
using Statistics
using GLM
using CairoMakie
using Random

# [Data Ingestion and Tier 1 Aggregation remain exactly the same as before...]

# ==========================================
# 3. Bootstrapped Cross-Validation Loop
# ==========================================
bin_sizes_to_test = 2:20
iterations = 100 # Number of bootstrap iterations per bin size
sample_fraction = 0.8 # Keep 80% of data per iteration

results = DataFrame(Bin_Size = Int[], Mean_True_SSR = Float64[], Std_True_SSR = Float64[])

println("Running Bootstrapped Cross-Validation Optimization (This may take a moment)...")

for w in bin_sizes_to_test
    iteration_errors = Float64[]
    
    for i in 1:iterations
        # Bootstrapping: Randomly sample 80% of the families
        df_family_boot = df_family[shuffle(1:nrow(df_family))[1:floor(Int, nrow(df_family)*sample_fraction)], :]
        
        # Bin the bootstrapped data
        df_family_boot.gap_bin = floor.(df_family_boot.family_gap_mean ./ w) .* w .+ (w / 2.0)
        
        df_binned = combine(groupby(df_family_boot, :gap_bin), 
            :log_family_norm_C => median => :log_overlap_median,
            nrow => :family_count
        )
        
        filter!(row -> row.family_count >= 2, df_binned)
        sort!(df_binned, :gap_bin)
        
        # Run regression
        horizon_bp, m_left, m_right = find_compaction_horizon(df_binned)
        
        if m_left === nothing || m_right === nothing
            continue # Skip failed alignments
        end
        
        # Calculate True Error against the BOOTSTRAPPED raw data
        df_left = filter(r -> r.family_gap_mean >= horizon_bp, df_family_boot)
        df_right = filter(r -> r.family_gap_mean < horizon_bp, df_family_boot)
        
        pred_left = predict(m_left, DataFrame(gap_bin = df_left.family_gap_mean))
        pred_right = predict(m_right, DataFrame(gap_bin = df_right.family_gap_mean))
        
        ssr_left = sum((df_left.log_family_norm_C .- pred_left).^2)
        ssr_right = sum((df_right.log_family_norm_C .- pred_right).^2)
        
        push!(iteration_errors, ssr_left + ssr_right)
    end
    
    # Calculate the mean and standard deviation of the error for this bin size
    if !isempty(iteration_errors)
        push!(results, (w, mean(iteration_errors), std(iteration_errors)))
    else
        push!(results, (w, Inf, 0.0))
    end
end

# Find the absolute best robust bin size
filter!(r -> r.Mean_True_SSR < Inf, results)
best_run = sort(results, :Mean_True_SSR)[1, :]

println("-------------------------------------------------")
println("🏆 ROBUST OPTIMAL BIN SIZE: ", best_run.Bin_Size, " bp")
println("-------------------------------------------------")

# ==========================================
# 4. Visualization of the Smoothed Curve
# ==========================================
fig = Figure(size = (800, 500), font = "Arial")
ax = Axis(fig[1, 1],
    title = "Bootstrapped Cross-Validation: Spatial Bin Size",
    subtitle = "Mean SSR over 50 iterations (80% subsampling) ± 1 Standard Deviation",
    xlabel = "Candidate Bin Width (bp)",
    ylabel = "Mean True Error (SSR)",
    titlesize = 18
)

# Plot an error band (Standard Deviation) to show stability
band!(ax, results.Bin_Size, 
      results.Mean_True_SSR .- results.Std_True_SSR, 
      results.Mean_True_SSR .+ results.Std_True_SSR, 
      color = (:dodgerblue, 0.2))

# Plot the smoothed mean error curve
lines!(ax, results.Bin_Size, results.Mean_True_SSR, color = :dodgerblue, linewidth = 3)
scatter!(ax, results.Bin_Size, results.Mean_True_SSR, color = :black, markersize = 8)

# Highlight the robust winner
scatter!(ax, [best_run.Bin_Size], [best_run.Mean_True_SSR], 
    color = :gold, markersize = 20, marker = :star5, strokewidth = 1, strokecolor = :black,
    label = "Robust Optimal Bin Size ($(best_run.Bin_Size) bp)")

axislegend(ax, position = :ct)
display(fig)
save("bootstrapped_bin_optimization.png", fig, px_per_unit = 2)

# ==========================================
# 6. Defining the Extreme Cutoff (3-Sigma Rule)
# ==========================================
# 1. Isolate the families in the Relaxed Baseline (>= 155 bp)
df_baseline = filter(r -> r.family_gap_mean >= 155.0, df_family)

# 2. Calculate the Mean and Standard Deviation of the baseline overlaps
baseline_mean = mean(df_baseline.log_family_norm_C)
baseline_std = std(df_baseline.log_family_norm_C)

# 3. Calculate the 3-Sigma Mathematical Ceiling
three_sigma_threshold = baseline_mean + (3.0 * baseline_std)

# 4. Find exactly where the Compaction Line crosses this ceiling
# The equation for the line is Y = mX + b. We solve for X: X = (Y - b) / m
m_compact = coef(model_compact)[2] # Slope
b_compact = coef(model_compact)[1] # Intercept

extreme_cutoff_bp = (three_sigma_threshold - b_compact) / m_compact

println("-------------------------------------------------")
println("📊 Baseline Mean (Log10): ", round(baseline_mean, digits=3))
println("📈 3-Sigma Threshold (Log10): ", round(three_sigma_threshold, digits=3))
println("🚨 EXTREME CRITICAL CUTOFF: ", round(extreme_cutoff_bp, digits=1), " bp")
println("-------------------------------------------------")

filter(row -> row.family_gap_mean <= 80, df_family)


