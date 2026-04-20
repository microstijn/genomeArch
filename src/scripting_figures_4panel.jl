using CairoMakie
using DataFrames
using CSV
using Statistics
using Colors # for robust color handling
using StatsBase

# ==============================================================================
# 1. LOAD AND PREPARE REAL DATA
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\genarch_with_full_taxonomy.csv", DataFrame)
println(names(df))
# Clean out critical missing data and tiny fragments
dropmissing!(df, [:absolute_gap_mean, :coding_density_pct, :mean_gene_length])
filter!(row -> !isnan(row.absolute_gap_mean) && !isnan(row.coding_density_pct), df)
df.total_genes = df.p_gene_nr .+ df.n_gene_nr
filter!(row -> row.total_genes > 100, df)

# Normalize count-based metrics to "per 1000 genes"
println("Calculating normalized metrics...")
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_boundary_collision = (df.boundary_collisions_nr ./ df.total_genes) .* 1000.0
df.norm_nested = (df.nested_genes_nr ./ df.total_genes) .* 1000.0
# Define head-to-head (Divergent) pairs specifically normalized
df.norm_Divergent_Pairs = (df.divergent_pairs_nr ./ df.total_genes) .* 1000.0

# Define a shared X-range for all plots
X_RANGE = (0.0, 500.0)
N_BINS = 80 # resolution of the raster bins

# ==============================================================================
# 2. CREATE THE LARGE RASTER MATRIX (4x3 Grid)
# ==============================================================================
println("Setting up Figure layout...")
# A very large canvas for publication quality
fig = Figure(size = (1800, 2200), fontsize = 18, font = "Helvetica")

# Common main title for the dashboard
Label(fig[0, 1:3], "The Macro-Biology of Genome Crush: Comprehensive Raster Analysis", fontsize = 30, font = :bold)
Label(fig[5, 1:3], "Mean Absolute Intergenic Space (bp) [Increasing Compaction →]", fontsize = 24, font = :bold)

# Use vibrant turbo colormap for maximum visibility of density differences
cmap = :turbo

# Define a function to easily add a hist2d panel with a log-scale colorbar
# Define a function to easily add a 2D histogram panel
using StatsBase # Ensure this is at the top of your script

function add_raster_panel!(pos, x_data, y_data, title_str, y_label, y_range)
    ax = Axis(pos, 
        title = title_str, 
        ylabel = y_label, 
        limits = (X_RANGE, y_range), 
        xreversed = true)
    
    # 1. Manually bin for the Heatmap (The Background Population)
    edges_x = range(X_RANGE[1], X_RANGE[2], length=N_BINS+1)
    edges_y = range(y_range[1], y_range[2], length=N_BINS+1)
    h_fit = fit(Histogram, (x_data, y_data), (edges_x, edges_y))
    
    # Plot the population density in a muted colormap so the frontier lines pop
    h = heatmap!(ax, h_fit.edges[1], h_fit.edges[2], log10.(h_fit.weights .+ 1); 
        colormap = :oslo) # 'oslo' is great for backgrounds; it's a subtle blue/black/white
    
    # 2. Calculate the Frontier (Moving Quantiles)
    # We want to see the 50th (median), 95th, and 99th percentiles per X-bin
    bin_centers = (edges_x[1:end-1] .+ edges_x[2:end]) ./ 2
    q50 = Float64[]
    q95 = Float64[]
    q99 = Float64[]
    
    for i in 1:(length(edges_x)-1)
        mask = edges_x[i] .<= x_data .< edges_x[i+1]
        if sum(mask) > 5
            vals = y_data[mask]
            push!(q50, median(vals))
            push!(q95, quantile(vals, 0.95))
            push!(q99, quantile(vals, 0.99))
        else
            # Fill with NaN if bin is too empty to keep line indices aligned
            push!(q50, NaN); push!(q95, NaN); push!(q99, NaN)
        end
    end

    # 3. Plot the Frontier Lines
    lines!(ax, bin_centers, q50, color = :white, linewidth = 1.5, label = "Median")
    lines!(ax, bin_centers, q95, color = :yellow, linewidth = 2, label = "95th Percentile")
    lines!(ax, bin_centers, q99, color = :red, linewidth = 2, label = "99th Percentile (The Frontier)")

    # 4. Reference Wall
    vlines!(ax, [134.0], color = (:white, 0.4), linestyle = :dash)
    
    
    return ax
end

println("Plotting RASTER rows...")
# ------------------------------------------------------------------------------
# ROW 1: Macro Global State
# ------------------------------------------------------------------------------
# Coding Density asymptotes to 100%
add_raster_panel!(fig[1, 1], df.absolute_gap_mean, df.coding_density_pct, 
    "R1C1. Global Coding Density", "Density (%)", (70.0, 100.0))

# Gene Length shows the "Protected" then "Truncated" phases
add_raster_panel!(fig[1, 2], df.absolute_gap_mean, df.mean_gene_length, 
    "R1C2. Gene Length (Shrinkage)", "Mean Gene Length (bp)", (500.0, 1400.0))

# Strand asymmetry shows if compaction forces genomes single-stranded (near 0.5 is balanced, 0.9+ is unbalanced)
add_raster_panel!(fig[1, 3], df.absolute_gap_mean, df.strand_asymmetry, 
    "R1C3. Strand Asymmetry", "Prop. (+) Genes", (0.3, 1.0))

# ------------------------------------------------------------------------------
# ROW 2: Operon Mechanics & Regulatory Space
# ------------------------------------------------------------------------------
# Operonicity increases linearly
add_raster_panel!(fig[2, 1], df.absolute_gap_mean, df.operonicity_score, 
    "R2C1. Operonicity Score", "% Genes in Operons", (0.0, 100.0))

# Mean Operon Size increases
add_raster_panel!(fig[2, 2], df.absolute_gap_mean, df.mean_operon_size, 
    "R2C2. Operon Cluster Size", "Mean Genes per Operon", (1.5, 6.0))

# True Inter-operon Regulatory Space (Should look like Absolute Gaps but shifted)
add_raster_panel!(fig[2, 3], df.absolute_gap_mean, df.inter_operon_gap_mean, 
    "R2C3. Inter-Operon Regulatory Space", "Mean regulatory gap (bp)", (0.0, 1500.0))

# ------------------------------------------------------------------------------
# ROW 3: Topological Response (The Spikes)
# ------------------------------------------------------------------------------
# Strand switches drop dramatically
add_raster_panel!(fig[3, 1], df.absolute_gap_mean, df.strand_switch_rate, 
    "R3C1. Strand Switch Frequency", "Switches per Gene", (0.0, 1.0))

# U-overlaps (same-strand) spike heavily
add_raster_panel!(fig[3, 2], df.absolute_gap_mean, df.norm_U, 
    "R3C2. Unidirectional Overlaps (U)", "U-Events per 1000 Genes", (0.0, 250.0))

# Boundary Collisions spike heavily
add_raster_panel!(fig[3, 3], df.absolute_gap_mean, df.norm_boundary_collision, 
    "R3C3. Boundary Collisions (C+D)", "Col-Events per 1000 Genes", (0.0, 300.0))

# ------------------------------------------------------------------------------
# ROW 4: Extreme and Avoided Topologies
# ------------------------------------------------------------------------------
# Divergent Head-to-Head Promoter collisions are *strictly* avoided (should stay low)
# Note: Divergent pairs metric in dashboard turn 16 was count, but divergent pairs is a specific Head-to-Head metric
# Let's normalize divergent pairs specifically
add_raster_panel!(fig[4, 1], df.absolute_gap_mean, df.norm_Divergent_Pairs, 
    "R4C1. Head-to-Head Pairs (D)", "D-Pairs per 1000 Genes", (0.0, 150.0))

# Entanglement Hubs explode vertically
add_raster_panel!(fig[4, 2], df.absolute_gap_mean, df.max_overlaps_per_gene, 
    "R4C2. Entanglement Hubs", "Max Overlaps per Gene", (0.0, 10.0))

# Nested Genes (total engulfment) are rare but appear only at extreme crush
add_raster_panel!(fig[4, 3], df.absolute_gap_mean, df.norm_nested, 
    "R4C3. Nested Genes", "Nest-Events per 1000 Genes", (0.0, 30.0))


# ==============================================================================
# 3. SAVE HIGH-RESOLUTION RESULT
# ==============================================================================
# Reduce the padding around colorbars to keep the grid tight

rowgap!(fig.layout, 10)
colgap!(fig.layout, 15)

output_path = raw"D:\pipeline_output\grand_compaction_raster_matrix.png"
save(output_path, fig, px_per_unit = 3) # px_per_unit=3 guarantees 300+ DPI result
println("Large raster matrix saved successfully to: $output_path")

display(fig)