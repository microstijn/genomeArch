# Setup Environment 
using Pkg
# Pkg.add(["CSV", "DataFrames", "CairoMakie", "Statistics", "GLM"]) # Uncomment if needed

using CSV
using DataFrames
using CairoMakie
using Statistics
using GLM # Required for piecewise regression in Panel C

println("Initializing publication figure generation script...")

# ==========================================================================
# 1. DATA LOADING & PREPROCESSING (Cumulative Data - Panels A, B)
# ==========================================================================
println("Loading data...")
# Load raw overlap pairs (exemplars and population)
df = CSV.File(raw"D:\pipeline_output\overlapping_pairs_annotated.csv") |> DataFrame
# Load full genome stats (gene counts, etc.)
dfs = CSV.read(raw"D:\pipeline_output\merged_imputed_ogt.csv", DataFrame)

opp_strand_df = filter(row -> row.overlap_type in ["Convergent", "Divergent"], df)

# --- ISOLATE EXEMPLARS VS. GENERAL POPULATION ---
println("Processing baseline distributions...")
exemplar_genomes = Set(filter(:overlap_length => >(120), opp_strand_df).genome)
exemplar_df = filter(:genome => g -> g in exemplar_genomes, opp_strand_df)

gene_counts = Dict(dfs.accession .=> dfs.p_gene_nr .+ dfs.n_gene_nr)
total_genes_pop = sum(get(gene_counts, g, 0) for g in unique(opp_strand_df.genome))
total_genes_ex = sum(get(gene_counts, g, 0) for g in exemplar_genomes)

# --- FAST CUMULATIVE CALCULATION & PERCENTAGE MISSING ---
max_x = 200 # Define max overlap length to plot
lengths = collect(1:max_x)

pop_counts = zeros(Int, max_x)
ex_counts = zeros(Int, max_x)

# Count overlaps by length
for x in 1:max_x
    pop_counts[x] = nrow(filter(:overlap_length => ==(x), opp_strand_df))
    ex_counts[x] = nrow(filter(:overlap_length => ==(x), exemplar_df))
end

# Calculate reverse cumulative sums
cum_pop_counts = reverse(cumsum(reverse(pop_counts)))
cum_ex_counts = reverse(cumsum(reverse(ex_counts)))

# Normalize (per 1000 genes)
cum_pop_density = (cum_pop_counts ./ total_genes_pop) .* 1000.0
cum_ex_density = (cum_ex_counts ./ total_genes_ex) .* 1000.0

# Log transform (log10(x + 1))
log_pop = log10.(cum_pop_density .+ 1.0)
log_ex  = log10.(cum_ex_density .+ 1.0)

# Calculate Percentage Missing
pct_missing = zeros(Float64, max_x)
for x in 1:max_x
    if cum_ex_density[x] > 0
        missing_val = ((cum_ex_density[x] - cum_pop_density[x]) / cum_ex_density[x]) * 100.0
        pct_missing[x] = max(0.0, missing_val) # Clamp to 0
    else
        pct_missing[x] = 0.0
    end
end

# ==========================================================================
# 2. DATA PREPROCESSING (Biophysical Data - Panel C - Example)
# ==========================================================================
# I have implemented the phylogenetically corrected piecewise logic here
# as a representative Panel C. Replace this section if another plot is intended.
println("Processing biophysical data for Panel C...")

# Prepare data for Panel C analysis
dfs.norm_C = (dfs.C_overlap_nr ./ (dfs.p_gene_nr .+ dfs.n_gene_nr)) .* 1000.0
dfc_clean = dropmissing(dfs, [:absolute_gap_mean, :norm_C, :family])
filter!(row -> 0 <= row.absolute_gap_mean <= 300, dfc_clean)

# Tier 1: Taxonomic Aggregation (collapse by family)
df_family = combine(groupby(dfc_clean, :family),
    :absolute_gap_mean => median => :family_gap_mean,
    :norm_C => median => :family_norm_C
)
df_family.log_family_norm_C = log10.(df_family.family_norm_C .+ 1.0)

# Tier 2: Spatial Binning (10-bp bins)
bin_width_c = 10
df_family.gap_bin = floor.(Int, df_family.family_gap_mean ./ bin_width_c) .* bin_width_c .+ (bin_width_c / 2.0)
df_binned = combine(groupby(df_family, :gap_bin),
    :log_family_norm_C => median => :log_overlap_median,
    nrow => :family_count
)
# Require at least 2 independent families per bin
filter!(row -> row.family_count >= 2, df_binned)
sort!(df_binned, :gap_bin)

# --- Define Piecewise Regression Function for Panel C ---
function find_c_horizon(df_bin)
    x = df_bin.gap_bin
    y = df_bin.log_overlap_median
    best_ssr, best_bp = Inf, 0
    best_m_left, best_m_right = nothing, nothing

    # Skip regression if too few data points
    if nrow(df_bin) < 8
        return best_bp, best_m_left, best_m_right
    end

    for i in 4:(length(x) - 4) # Avoid testing too close to the edges
        bp = x[i]
        df_left = filter(r -> r.gap_bin >= bp, df_bin)
        df_right = filter(r -> r.gap_bin < bp, df_bin)

        # Skip if a split doesn't have enough points for a line
        if nrow(df_left) < 3 || nrow(df_right) < 3
            continue
        end

        m_left = lm(@formula(log_overlap_median ~ gap_bin), df_left)
        m_right = lm(@formula(log_overlap_median ~ gap_bin), df_right)
        total_ssr = deviance(m_left) + deviance(m_right)

        if total_ssr < best_ssr
            best_ssr = total_ssr
            best_bp = bp
            best_m_left = m_left
            best_m_right = m_right
        end
    end
    return best_bp, best_m_left, best_m_right
end

# Calculate the optimal horizon breakpoint
println("Calculating Panel C horizon breakpoint...")
horizon_bp_c, model_left_c, model_right_c = find_c_horizon(df_binned)
println("-------------------------------------------------")
println("🚀 PANEL C CALCULATED COMPACTION HORIZON: ", horizon_bp_c, " bp")
println("-------------------------------------------------")

# ==========================================================================
# 3. VISUALIZATION (PAPER-INSPIRED THEME & LAYOUT)
# ==========================================================================
println("Generating plot...")

# Define the minimalist, paper-inspired theme
paper_theme = Theme(
    font = "Helvetica",
    fontsize = 14,
    Axis = (
       xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,      #Open top
        rightspinevisible = false,    #Open right
        spinewidth = 1.5,            #Thicker, bolder axes
        bottomspinecolor = :black,
        leftspinecolor = :black,
        xtickwidth = 1.5,
        ytickwidth = 1.5,
    )
)
set_theme!(paper_theme)

# Reordered layout definition: Figure size (1200, 450)
fig = Figure(size = (1200, 450), backgroundcolor = :white)

# Define a central colorbar axis (shared color scale for hexbin/medians)
# This is currently defined inside Panel C, but a shared one could be placed in fig[1,4] if needed.

# --- PANEL A: Log-Cumulative Distribution ---
# Moved to fig[1,1] - Add xlabel as it's no longer shared in a vertical stack
ax1 = Axis(fig[1, 1],
    xlabel = "Overlap length (bp)",
    ylabel = "Cumulative overlaps\n(log₁₀ per 1000 genes)"
)

Label(fig[1, 1, TopLeft()], "a", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

vlines!(ax1, [120], color = :black, linestyle = :dash, linewidth = 1.5)

# Shaded Area
band_x = 1:max_x
band_lower = log_pop
band_upper = max.(log_pop, log_ex)
band!(ax1, band_x, band_lower, band_upper, color = (:grey70, 0.5), label = "missing overlap space")

# Lines
lines!(ax1, lengths, log_pop, color = :grey60, linewidth = 2.5, label = "all genomes (filtered)")
lines!(ax1, lengths, log_ex, color = :mediumseagreen, linewidth = 3.0, label = "exemplars (unfiltered)")

ylims!(ax1, 0, maximum(log_ex) * 1.1)

# Clean legend, no border box
axislegend(ax1, position = :rt, framevisible = false, labelsize=12)


# --- PANEL B: Percentage Missing Gradient ---
# Moved to fig[1,2]
ax2 = Axis(fig[1, 2],
xlabel = "Overlap length (bp)", # Xlabel is now independent
ylabel = "Estimated missing (%)"
)

Label(fig[1, 2, TopLeft()], "b", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

vlines!(ax2, [120], color = :black, linestyle = :dash, linewidth = 1.5)

band!(ax2, lengths, zeros(max_x), pct_missing, color = (:mediumpurple, 0.3))
lines!(ax2, lengths, pct_missing, color = :mediumpurple, linewidth = 3.0)

ylims!(ax2, 0, 105) 


# --- PANEL C: Biophysical Compression (Piecewise Regression) ---
# Moved to fig[1,3]
# **NOTE TO USER: PLACE YOUR PANEL C PLOTTING CODE HERE IF IT DIFFERS**

# Define Axis 3
ax3 = Axis(fig[1, 3],
    xlabel = "Abs. gap mean (bp) [Decreasing →]",
    ylabel = "Log₁₀(Convergent overlaps + 1)",
    xreversed = true # Decreasing gap size to the right
)
Label(fig[1, 3, TopLeft()], "c", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

# I am plotting the robust Hexbin + Medians view discussed previously
# 1. Hexbin plotting the underlying Family Data
hexbin!(ax3, df_family.family_gap_mean, df_family.log_family_norm_C,
    cellsize = (10, 0.1), 
    colormap = cgrad(:Blues, rev=false),
    strokewidth = 0.5,
    strokecolor = :white
)

# 2. Plot the 10-bp Binned Medians (The data driving the math)
scatter!(ax3, df_binned.gap_bin, df_binned.log_overlap_median, 
    color = :darkorange, markersize = 12, label = "10-bp Binned Medians (Families)", 
    strokecolor = :black, strokewidth = 1)

# 3. Plot the piecewise regression lines
x_relaxed = Float64.(filter(x -> x >= horizon_bp_c, df_binned.gap_bin))
x_compact = Float64.(filter(x -> x < horizon_bp_c, df_binned.gap_bin))
y_relaxed_pred = Float64.(predict(model_left_c, DataFrame(gap_bin = x_relaxed)))
y_compact_pred = Float64.(predict(model_right_c, DataFrame(gap_bin = x_compact)))

lines!(ax3, x_relaxed, y_relaxed_pred, color = :red, linewidth = 4, label = "Relaxed Fit")
lines!(ax3, x_compact, y_compact_pred, color = :darkred, linewidth = 4, label = "Compaction Fit")

# 4. Add the vertical Horizon boundary (labeled in label, not on plot)
vlines!(ax3, [horizon_bp_c], color = :black, linestyle = :dash, linewidth = 2, 
    label = "Compaction Horizon ($horizon_bp_c bp)")

# Formatting
xlims!(ax3, 400, 0)
# Use top padding for limit rather than absolute number to match A/B
ylims!(ax3, -0.05, maximum(df_family.log_family_norm_C) * 1.1) 
axislegend(ax3, position = :lt, framevisible=false, labelsize=12) # Top-Left legend matching other panels

# --- FINAL LAYOUT FORMATTING ---
# linkxaxes!(ax1, ax2) # CANNOT LINK AXES SIDE-BY-SIDE
xlims!(ax1, 0, max_x)
xlims!(ax2, 0, max_x)
xlims!(ax3, 300, 0)
# Adjust spacing between columns
colgap!(fig.layout, 25)

display(fig)
save("combined_overlap_publication_figure.png", fig, px_per_unit = 4) # Higher resolution export
save("combined_overlap_publication_figure.pdf", fig)
println("Done! Publication-styled reordered figure saved.")