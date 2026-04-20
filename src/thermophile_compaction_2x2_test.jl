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
# Calculate genome size in Megabases for cleaner X-axes
df.genome_size_mb = df.genome_size ./ 1_000_000.0

dropmissing!(df, [:absolute_gap_mean, :genome_size, :OGT])
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
# Increased height to 900 to comfortably fit two rows
fig = Figure(size = (1200, 900), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "Thermophile Compaction: Overlaps Driven by Physical Space Limits", fontsize = 22, font = :bold)

axes = []

# Define the 2x2 grid structure: (Row, Col, X-var, Y-var, Title, X-label, Y-label, Color)
panels = [
    # ROW 1: Gap Mean
    (1, 1, :absolute_gap_mean, :norm_U, "1. Unidirectional Overlaps", "Absolute Gap Mean (bp) [Decreasing →]", "Events / 1k Genes", :steelblue),
    (1, 2, :absolute_gap_mean, :norm_C, "2. Convergent Overlaps", "Absolute Gap Mean (bp) [Decreasing →]", "Events / 1k Genes", :firebrick),
    
    # ROW 2: Genome Size
    (2, 1, :genome_size_mb, :norm_U, "3. Unidirectional Overlaps", "Genome Size (Mb) [Decreasing →]", "Events / 1k Genes", :steelblue),
    (2, 2, :genome_size_mb, :norm_C, "4. Convergent Overlaps", "Genome Size (Mb) [Decreasing →]", "Events / 1k Genes", :firebrick)
]

# ==============================================================================
# 3. PLOT THE 2x2 GRID
# ==============================================================================
for (row, col, x_sym, y_sym, title_str, x_label, y_label, color) in panels
    # Force X-limits for Gap Mean to match your original script
    x_lims = x_sym == :absolute_gap_mean ? (0, 400) : (nothing, nothing)
    
    ax = Axis(fig[row, col], 
              title = title_str, titlealign = :left,
              xlabel = x_label,
              ylabel = y_label,
              xreversed = true, # Time/Compaction flows left to right
              limits = (x_lims, (nothing, nothing)))
    push!(axes, ax)
    
    x_data = Float64.(df_thermo[!, x_sym])
    y_data = Float64.(df_thermo[!, y_sym])
    
    # Scatter plot of the thermophilic genomes
    scatter!(ax, x_data, y_data, color = (color, 0.4), markersize = 8, strokewidth = 0.5, strokecolor = :black)
    
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

# Link axes intelligently (link Y-axes for columns to compare overlaps directly)
linkyaxes!(axes[1], axes[3]) # Link Unidirectional Y-axes
linkyaxes!(axes[2], axes[4]) # Link Convergent Y-axes

colgap!(fig.layout, 20)
rowgap!(fig.layout, 20)

# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\thermophile_compaction_2x2_test.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)