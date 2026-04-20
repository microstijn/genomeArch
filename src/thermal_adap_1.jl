using CairoMakie
using DataFrames
using CSV
using Loess

# Load data
df = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

# Clean missing and NaN values
dropmissing!(df, [:absolute_gap_mean, :coding_density_pct, :OGT])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.OGT), df)

# Setup figure
fig = Figure(size = (1200, 500), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Thermal Adaptation: Nucleoid Volume Compression", fontsize = 24, font = :bold)

# Panel 1: Physical Compaction Wall (Absolute Gaps)
ax1 = Axis(fig[1, 1], 
           title = "1. Intergenic Space vs Temperature", 
           xlabel = "Optimal Growth Temperature (°C)", 
           ylabel = "Absolute Gap Mean (bp)")

scatter!(ax1, df.OGT, df.absolute_gap_mean, color = (:firebrick, 0.4), markersize = 8)

# Panel 2: Total Coding Density Limit
ax2 = Axis(fig[1, 2], 
           title = "2. Coding Density vs Temperature", 
           xlabel = "Optimal Growth Temperature (°C)", 
           ylabel = "Coding Density (%)")

scatter!(ax2, df.OGT, df.coding_density_pct, color = (:steelblue, 0.4), markersize = 8)

# Add smoothing lines to detect phase transitions across temperatures
for (ax, y_col) in zip([ax1, ax2], [:absolute_gap_mean, :coding_density_pct])
    model = loess(df.OGT, df[!, y_col], span=0.5)
    x_eval = range(minimum(df.OGT), maximum(df.OGT), length=100)
    y_eval = predict(model, x_eval)
    lines!(ax, x_eval, y_eval, linewidth=4, color=:black)
end

save(raw"D:\pipeline_output\thermal_adaptation_ogt_only.png", fig, px_per_unit = 3)
display(fig)