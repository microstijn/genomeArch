using CairoMakie
using DataFrames
using CSV
using StatsBase
using LinearAlgebra
using Random

# ==============================================================================
# 1. THE DISTANCE CORRELATION (dCor) ALGORITHM
# ==============================================================================
function distance_correlation(x::AbstractVector, y::AbstractVector)
    n = length(x)
    a = abs.(x .- x')
    b = abs.(y .- y')
    
    a_row_mean, a_col_mean, a_mean = mean(a, dims=2), mean(a, dims=1), mean(a)
    A = a .- a_row_mean .- a_col_mean .+ a_mean
    
    b_row_mean, b_col_mean, b_mean = mean(b, dims=2), mean(b, dims=1), mean(b)
    B = b .- b_row_mean .- b_col_mean .+ b_mean
    
    dcov_xy = sqrt(sum(A .* B) / n^2)
    dvar_x = sqrt(sum(A .* A) / n^2)
    dvar_y = sqrt(sum(B .* B) / n^2)
    
    return (dvar_x * dvar_y == 0) ? 0.0 : dcov_xy / sqrt(dvar_x * dvar_y)
end

# ==============================================================================
# 2. LOAD AND PREPARE NUMERIC MATRIX (Counts & Lengths)
# ==============================================================================
println("Loading dataset...")
df = CSV.read(raw"D:\pipeline_output\genarch_with_full_taxonomy.csv", DataFrame)

# Calculate normalized overlap counts
df.total_genes = df.p_gene_nr .+ df.n_gene_nr
df.norm_U = ((df.p_U_overlap_nr .+ df.n_U_overlap_nr) ./ df.total_genes) .* 1000.0
df.norm_C = (df.C_overlap_nr ./ df.total_genes) .* 1000.0
df.norm_D = (df.D_overlap_nr ./ df.total_genes) .* 1000.0
df.norm_abutting = (df.abutting_genes_nr ./ df.total_genes) .* 1000.0
df.norm_nested = (df.nested_genes_nr ./ df.total_genes) .* 1000.0

# Calculate Mean Overlap Lengths
# We use max.(1.0, count) in the denominator to prevent DivideByZero (NaN) when count is 0
df.mean_U_length = (df.p_U_overlap_length_sum .+ df.n_U_overlap_length_sum) ./ max.(1.0, df.p_U_overlap_nr .+ df.n_U_overlap_nr)
df.mean_C_length = df.C_length_sum ./ max.(1.0, df.C_overlap_nr)
df.mean_D_length = df.D_length_sum ./ max.(1.0, df.D_overlap_nr)

# Zero out the length if the count was actually 0
df.mean_U_length[df.p_U_overlap_nr .+ df.n_U_overlap_nr .== 0] .= 0.0
df.mean_C_length[df.C_overlap_nr .== 0] .= 0.0
df.mean_D_length[df.D_overlap_nr .== 0] .= 0.0

# The strictly non-circular topological variables + OVERLAP LENGTHS
vars = [
    :absolute_gap_mean, :coding_density_pct, :mean_gene_length, 
    :strand_switch_rate, 
    :norm_U, :norm_C, :norm_D, 
    :mean_U_length, :mean_C_length, :mean_D_length,
    :norm_abutting, :norm_nested, :max_overlaps_per_gene
]

df_clean = dropmissing(df[:, vars])
filter!(row -> all(!isnan, row), df_clean)
n_vars = length(vars)
labels = string.(vars)

# ==============================================================================
# 3. STRATIFIED SUBSAMPLING & MATRIX CALCULATION
# ==============================================================================
println("Performing Stratified Subsampling along the Compaction Gradient...")
Random.seed!(42) 

df_clean.original_idx = 1:nrow(df_clean)
sampled_indices = Int[]

bin_edges = 0.0:15.0:500.0
MAX_PER_BIN = 60 

for i in 1:(length(bin_edges)-1)
    mask = bin_edges[i] .<= df_clean.absolute_gap_mean .< bin_edges[i+1]
    indices_in_bin = df_clean.original_idx[mask]
    n_avail = length(indices_in_bin)
    
    if n_avail > 0
        n_take = min(n_avail, MAX_PER_BIN)
        append!(sampled_indices, sample(indices_in_bin, n_take, replace=false))
    end
end

println("Total genomes selected for dCor calculation: ", length(sampled_indices))

raw_data_matrix = Matrix(Float64.(df_clean[sampled_indices, 1:n_vars]))

valid_cols = [var(raw_data_matrix[:, i]) > 0.0 for i in 1:n_vars]
data_matrix = raw_data_matrix[:, valid_cols]
labels = labels[valid_cols]
n_vars = length(labels) 

if sum(.!valid_cols) > 0
    println("Dropped flatline variables: ", vars[.!valid_cols])
end

println("Calculating Spearman Matrix...")
spearman_mat = corspearman(data_matrix)

println("Calculating Distance Correlation (dCor) Matrix...")
dcor_mat = zeros(n_vars, n_vars)

for i in 1:n_vars
    for j in 1:n_vars
        if i == j
            dcor_mat[i, j] = 1.0
        else
            if j > i 
                val = distance_correlation(data_matrix[:, i], data_matrix[:, j])
                dcor_mat[i, j] = val
                dcor_mat[j, i] = val
            end
        end
    end
end

println("Calculating Delta Matrix...")
delta_mat = dcor_mat .- abs.(spearman_mat)

# ==============================================================================
# 4. PLOT THE 3-PANEL DASHBOARD
# ==============================================================================
fig = Figure(size = (2000, 700), fontsize = 14, font = "Helvetica")

function add_heatmap!(pos, matrix, title_str, cmap, color_range)
    ax = Axis(pos, title = title_str, 
        xticks = (1:n_vars, labels), yticks = (1:n_vars, labels),
        xticklabelrotation = pi/4)
    
    h = heatmap!(ax, 1:n_vars, 1:n_vars, matrix', 
        colormap = cmap, colorrange = color_range)
    
    Colorbar(pos[1, 2], h)
    return ax
end

add_heatmap!(fig[1, 1], abs.(spearman_mat), "A. Absolute Spearman |ρ|", :Blues, (0.0, 1.0))
add_heatmap!(fig[1, 2], dcor_mat, "B. Distance Correlation (dCor)", :Purples, (0.0, 1.0))
add_heatmap!(fig[1, 3], delta_mat, "C. The Delta (dCor - |ρ|)", :inferno, (0.0, maximum(delta_mat)))

output_path = raw"D:\pipeline_output\strict_topology_with_lengths_dcor.png"
save(output_path, fig, px_per_unit=3)
println("Dashboard saved to: $output_path")
display(fig)