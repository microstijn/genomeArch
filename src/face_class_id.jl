using CairoMakie
using DataFrames
using CSV

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
println("Loading dataset with exact TaxIDs...")
df = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

df.total_genes = df.p_gene_nr .+ df.n_gene_nr
df.genome_size_mb = df.genome_size ./ 1_000_000.0

dropmissing!(df, [:absolute_gap_mean, :genome_size, :class])
filter!(row -> !isnan(row.absolute_gap_mean) && row.total_genes > 100, df)

# Using your EXACT taxonomy IDs and labels
target_classes = [
    ("1236", "Gammaproteobacteria", :goldenrod),
    ("28211", "Alphaproteobacteria", :darkorange),
    ("91061", "Bacilli", :seagreen),
    ("1760", "Actinomycetes", :indigo),
    ("31969", "Mollicutes", :dodgerblue),
    ("204429", "Chlamydiia", :purple),
    ("188708", "Thermotogae", :firebrick),
    ("3028117", "Cyanophyceae", :teal)  
]

# ==============================================================================
# 2. SETUP THE FACETED FIGURE
# ==============================================================================
fig = Figure(size = (1800, 900), fontsize = 16, font = "Helvetica")
Label(fig[0, 1:4], "Evolutionary Trajectories per Class: Decay vs. Compaction", fontsize = 26, font = :bold)

# Loop through the classes and create a 2x4 grid of subplots
for (i, (tax_id, cls_name, cls_color)) in enumerate(target_classes)
    row = div(i - 1, 4) + 1
    col = mod1(i, 4)
    
    ax = Axis(fig[row, col], 
              title = cls_name, titlealign = :center,
              xlabel = "Absolute Gap Mean (bp) [Decreasing →]",
              ylabel = col == 1 ? "Genome Size (Mb)" : "", 
              xreversed = true,
              limits = ((0, 400), (0, 10)))
    
    # -- 1. Plot Background (All Data) --
    scatter!(ax, df.absolute_gap_mean, df.genome_size_mb, 
             color = (:grey, 0.05), markersize = 4, strokewidth = 0)
    
    # -- 2. Plot Target Class by EXACT TaxID --
    id_regex = Regex("^$(tax_id)\\b")
    sub_df = filter(row -> !ismissing(row.class) && occursin(id_regex, row.class), df)
    
    if nrow(sub_df) > 0
        scatter!(ax, sub_df.absolute_gap_mean, sub_df.genome_size_mb, 
                 color = (cls_color, 0.6), markersize = 6, strokewidth = 0.5, strokecolor = :black)
    end
end

# ==============================================================================
# 3. SAVE AND DISPLAY
# ==============================================================================
colgap!(fig.layout, 15)
rowgap!(fig.layout, 15)

output_path = raw"D:\pipeline_output\faceted_class_exact_ids.png"
save(output_path, fig, px_per_unit = 3)
println("Plot saved to: $output_path")
display(fig)