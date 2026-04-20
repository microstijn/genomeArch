using CairoMakie
using DataFrames
using CSV
using Loess
using Statistics

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset for Macro-Deletions...")
df = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

df.total_genes = df.p_gene_nr .+ df.n_gene_nr

# Ensure required columns exist and are valid
dropmissing!(df, [:absolute_gap_mean, :OGT, :genome_size, :total_genes])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.OGT) && row.total_genes > 100, df)

# Convert genome size to Megabases (Mb) for readability
df.genome_size_mb = df.genome_size ./ 1_000_000.0

# Define the Temperature Cohorts
cohorts = [
    ("Mesophiles (< 30°C)", row -> row.OGT < 30.0, :steelblue),
    ("Thermophiles (≥ 45°C)", row -> row.OGT >= 45.0, :firebrick)
]

# ==============================================================================
# 2. SETUP THE FIGURE
# ==============================================================================
fig = Figure(size = (1200, 500), fontsize = 18, font = "Helvetica")
Label(fig[0, 1:2], "The True Compaction Ledger: Mass Gene Deletion", fontsize = 24, font = :bold)

panels = [
    (:genome_size_mb, "1. Total Genome Size Collapse", "Genome Size (Mb)", 1),
    (:total_genes, "2. Total Gene Count Collapse", "Number of Genes", 2)
]

axes = []
eval_points = collect(400.0:-1.0:0.0)
bin_edges = 0.0:15.0:400.0
bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2.0

# ==============================================================================
# 3. PLOT COHORT TRAJECTORIES
# ==============================================================================
for (sym, title_str, y_label, col) in panels
    ax = Axis(fig[1, col], 
              title = title_str, titlealign = :left,
              xlabel = "Absolute Gap Mean (bp) [Decreasing →]",
              ylabel = y_label,
              xreversed = true, 
              limits = ((0, 400), (nothing, nothing)))
    push!(axes, ax)
    
    for (label_str, filter_func, color) in cohorts
        sub_df = filter(filter_func, df)
        
        x_data = sub_df.absolute_gap_mean
        y_data = Float64.(sub_df[!, sym])
        
        valid_centers, bin_means = Float64[], Float64[]
        for i in 1:(length(bin_edges)-1)
            mask = bin_edges[i] .<= x_data .< bin_edges[i+1]
            if sum(mask) >= 5 # Dynamic Bounding
                push!(valid_centers, bin_centers[i])
                push!(bin_means, mean(y_data[mask]))
            end
        end
        
        if length(valid_centers) > 3
            model = loess(valid_centers, bin_means, span=0.45)
            safe_eval = filter(x -> minimum(valid_centers) <= x <= maximum(valid_centers), eval_points)
            y_smooth = predict(model, safe_eval)
            
            # Plot the smoothed trends
            lines!(ax, safe_eval, y_smooth, linewidth = 5, color = color, label = label_str)
        end
    end
    
    if col == 1
        axislegend(ax, position = :lt, framevisible = false, labelsize = 16)
    end
end

linkxaxes!(axes...)
colgap!(fig.layout, 20)

# ==============================================================================
# 4. SAVE AND DISPLAY
# ==============================================================================
output_path = raw"D:\pipeline_output\macro_deletion_ledger.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)