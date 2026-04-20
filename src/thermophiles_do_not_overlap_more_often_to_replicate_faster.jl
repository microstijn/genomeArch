using CairoMakie
using DataFrames
using CSV
using Loess

# ==============================================================================
# 1. LOAD DATA AND STRIP IMPUTATIONS
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\master_imputed_genarch_with_madin.csv", DataFrame)

# Calculate derived metrics
df.total_genes = df.p_gene_nr .+ df.n_gene_nr
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.genome_size_mb = df.genome_size ./ 1_000_000.0

# THE CRITICAL STEP: Keep ONLY the genomes with direct laboratory measurements
filter!(row -> !ismissing(row.doubling_h_imputation_level) && 
               row.doubling_h_imputation_level == "direct", df)

println("Stripped all imputations. Retained $(nrow(df)) strictly empirical genomes.")

# Clean missing values for the plot
dropmissing!(df, [:norm_U, :doubling_h_imputed, :OGT])

filter!(:doubling_h_imputed => x -> x.> 0.0, df)

# Split into Thermophiles and Non-Thermophiles for clear biological comparison
df_thermo = filter(row -> row.OGT >= 45.0, df)
df_meso   = filter(row -> row.OGT < 45.0, df)

# ==============================================================================
# 2. SETUP THE FIGURE
# ==============================================================================
fig = Figure(size = (1200, 600), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Strict Empirical Test: Overlaps vs Speed (No Imputed Data)", fontsize = 22, font = :bold)

# Panel 1: Mesophiles (The Control)
ax1 = Axis(fig[1, 1], 
           title = "Non-Thermophiles (< 45°C)",
           xlabel = "Unidirectional Overlaps (per 1k genes)", 
           ylabel = "Doubling Time (Hours) [Faster ↑]")

# Panel 2: Thermophiles (The Wright et al. Test Group)
ax2 = Axis(fig[1, 2], 
           title = "Thermophiles (≥ 45°C)",
           xlabel = "Unidirectional Overlaps (per 1k genes)")

# ==============================================================================
# 3. PLOT EMPIRICAL DATA
# ==============================================================================
# Scatter the non-thermophiles (Blue)
scatter!(ax1, df_meso.norm_U, df_meso.doubling_h_imputed, 
         color = (:steelblue, 0.6), markersize = 8, strokewidth = 0.5, strokecolor = :black)

# Scatter the thermophiles (Red)
scatter!(ax2, df_thermo.norm_U, df_thermo.doubling_h_imputed, 
         color = (:firebrick, 0.8), markersize = 8, strokewidth = 0.5, strokecolor = :black)


linkyaxes!(ax1, ax2)
colgap!(fig.layout, 15)

for ax in (ax1, ax2)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.yscale = log10

end



# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\strict_empirical_growth_vs_overlaps.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)