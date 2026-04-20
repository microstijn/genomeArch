using CairoMakie
using DataFrames
using CSV

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\master_imputed_genarch_with_madin.csv", DataFrame)

# Calculate derived metrics
df.total_genes = df.p_gene_nr .+ df.n_gene_nr
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_C = (df.C_overlap_nr ./ df.total_genes) .* 1000.0

# Clean topological and thermal NaNs globally
dropmissing!(df, [:absolute_gap_mean, :norm_U, :norm_C, :OGT])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.OGT), df)

# ==============================================================================
# 2. SPLIT FOREGROUND (EMPIRICAL) AND BACKGROUND (UNKNOWN)
# ==============================================================================
# Foreground: Direct doubling time matches
df_empirical = filter(row -> !ismissing(row.doubling_h_imputation_level) && 
                             row.doubling_h_imputation_level == "direct" &&
                             !ismissing(row.doubling_h_imputed) && 
                             row.doubling_h_imputed > 0, df)

# Background: Everything else
df_background = filter(row -> ismissing(row.doubling_h_imputation_level) || 
                              row.doubling_h_imputation_level != "direct", df)

println("Foreground (Empirical Matches): $(nrow(df_empirical)) genomes.")
println("Background (Grey Universe): $(nrow(df_background)) genomes.")

# Split by temperature
df_bg_meso = filter(row -> row.OGT < 45.0, df_background)
df_bg_thermo = filter(row -> row.OGT >= 45.0, df_background)

df_emp_meso = filter(row -> row.OGT < 45.0, df_empirical)
df_emp_thermo = filter(row -> row.OGT >= 45.0, df_empirical)

# ==============================================================================
# 3. BIN THE BACKGROUND DATA FOR RAINCLOUDS
# ==============================================================================
# We define a bin width of 20 base pairs. We cast them to Ints so Makie's 
# rainclouds! treats them as categories but plots them at the correct X coordinates.
bin_width = 20.0
df_bg_meso.gap_cat = floor.(Int, df_bg_meso.absolute_gap_mean ./ bin_width) .* Int(bin_width) .+ Int(bin_width / 2)
df_bg_thermo.gap_cat = floor.(Int, df_bg_thermo.absolute_gap_mean ./ bin_width) .* Int(bin_width) .+ Int(bin_width / 2)

# ==============================================================================
# 4. SETUP THE FIGURE AND COLORMAP
# ==============================================================================
fig = Figure(size = (1600, 1000), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Architectural Rainclouds Colored by Empirical Speed (Log10)", fontsize = 24, font = :bold)

# Log10 color constraints for the empirical data
raw_doubling = Float64.(df_empirical.doubling_h_imputed)
c_min = minimum(raw_doubling)
c_max = quantile(raw_doubling, 0.95)
cmap = Reverse(:inferno) 

# --- ROW 1: Unidirectional Overlaps ---
ax11 = Axis(fig[1, 1], title = "Non-Thermophiles (< 45°C)",
            ylabel = "Unidirectional Overlaps\n(per 1k genes)",
            xreversed = true, limits = ((0, 400), (0, 400)))

ax12 = Axis(fig[1, 2], title = "Thermophiles (≥ 45°C)",
            xreversed = true, limits = ((0, 400), (0, 400       )))

# --- ROW 2: Convergent Overlaps ---
ax21 = Axis(fig[2, 1], xlabel = "Absolute Gap Mean (bp) [Decreasing →]", 
            ylabel = "Convergent Overlaps\n(per 1k genes)",
            xreversed = true, limits = ((0, 400), (0, 200)))

ax22 = Axis(fig[2, 2], xlabel = "Absolute Gap Mean (bp) [Decreasing →]", 
            xreversed = true, limits = ((0, 400), (0, 200   )))

# ==============================================================================
# 5. PLOT THE DATA
# ==============================================================================
function plot_native_raincloud!(ax, bg_data, emp_data, y_sym)
    
    # 1. The Background Universe (Using native rainclouds!)
    # This automatically handles the cloud, the boxplot, and the faint background rain.
    rainclouds!(ax,
                bg_data.gap_cat,
                Float64.(bg_data[!, y_sym]), 
                plot_boxplots = true, 
                boxplot_width = 5,
                center_boxplot = false,
                cloud_width = 20,
                #violin_limits = extrema,
                #gap = 0.5,
                hist_bins = 50,
                clouds= hist,
                #color = (:lightgrey, 0.4),   # Cloud/Box color
                markersize = 3,              # Faint jittered rain for the background
                side = :left) 

    # 2. The Empirical Foreground
    # We scatter the known speeds at their exact continuous coordinates so they float over the rainclouds.
    sc = scatter!(ax, emp_data.absolute_gap_mean, Float64.(emp_data[!, y_sym]), 
                  color = emp_data.doubling_h_imputed, 
                  colormap = cmap, colorrange = (c_min, c_max), colorscale = log10,
                  markersize = 10, strokewidth = 0.5, strokecolor = :white)
    return sc
end

# Draw the 4 panels
plot_native_raincloud!(ax11, df_bg_meso, df_emp_meso, :norm_U)
sc_ref = plot_native_raincloud!(ax12, df_bg_thermo, df_emp_thermo, :norm_U)

plot_native_raincloud!(ax21, df_bg_meso, df_emp_meso, :norm_C)
plot_native_raincloud!(ax22, df_bg_thermo, df_emp_thermo, :norm_C)

# Link axes strategically
linkyaxes!(ax11, ax12) 
linkyaxes!(ax21, ax22) 
linkxaxes!(ax11, ax12, ax21, ax22)

colgap!(fig.layout, 15)
rowgap!(fig.layout, 15)

# ==============================================================================
# 6. ADD LOG10 COLORBAR
# ==============================================================================
Colorbar(fig[1:2, 3], sc_ref, 
         label = "Empirical Doubling Time (Hours)\n[Log10 Scale | Brighter = Faster Growth]", 
         ticklabelsize = 16)

# ==============================================================================
# 7. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\topological_landscape_2x2_native_rainclouds.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)