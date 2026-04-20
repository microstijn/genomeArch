using CairoMakie
using DataFrames
using CSV
using Loess

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

df.total_genes = df.p_gene_nr .+ df.n_gene_nr
dropmissing!(df, [:absolute_gap_mean, :OGT])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.OGT) && row.total_genes > 100, df)

# Normalize topological counts per 1000 genes
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_C = (df.C_overlap_nr ./ df.total_genes) .* 1000.0

# Isolate the high-temperature population (e.g., OGT >= 45°C)
df_thermo = filter(row -> row.OGT >= 45.0, df)
println("Isolated $(nrow(df_thermo)) thermophilic genomes for analysis.")

# ==============================================================================
# 2. SETUP THE FIGURE
# ==============================================================================
fig = Figure(size = (1200, 500), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Thermophile Compaction: Overlaps Dictated by Physical Limits, Not Temperature", fontsize = 22, font = :bold)

axes = []
panels = [
    (:norm_U, "1. Unidirectional Overlaps (OGT ≥ 45°C)", "Events / 1k Genes", 1, :steelblue),
    (:norm_C, "2. Convergent Overlaps (OGT ≥ 45°C)", "Events / 1k Genes", 2, :firebrick)
]

# ==============================================================================
# 3. PLOT THERMOPHILES VS EVOLUTIONARY TIMELINE
# ==============================================================================
for (sym, title_str, y_label, col, color) in panels
    ax = Axis(fig[1, col], 
              title = title_str, titlealign = :left,
              xlabel = "Absolute Gap Mean (bp) [Decreasing →]",
              ylabel = y_label,
              xreversed = true, # Time flows left to right
              limits = ((0, 400), (nothing, nothing)))
    push!(axes, ax)
    
    x_data = df_thermo.absolute_gap_mean
    y_data = Float64.(df_thermo[!, sym])
    
    # Scatter plot of the thermophilic genomes
    scatter!(ax, x_data, y_data, color = (color, 0.4), markersize = 10)
    
    # Fit a standard Loess to visualize the phase transition
    try
        model = loess(x_data, y_data, span=0.5)
        x_eval = range(minimum(x_data), maximum(x_data), length=100)
        y_smooth = predict(model, x_eval)
        lines!(ax, x_eval, y_smooth, linewidth = 4, color = :black)
    catch e
        # Pass if data is too sparse for loess
    end
end

linkxaxes!(axes...)
colgap!(fig.layout, 20)

# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\thermophile_compaction_test.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)