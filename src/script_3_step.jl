using DataFrames
using GLM
using Distributions
using Statistics
using GaussianMixtures

println("======================================================")
println("=== FINAL UNIFIED MACRO-EVOLUTIONARY PIPELINE =======")
println("======================================================")

# ---------------------------------------------------------
# SETUP: Clean the data
# ---------------------------------------------------------
df_final = dropmissing(df_s, [:mean_gap_size, :C_overlap_density, :genome_size, :total_genes])
filter!(row -> !isnan(row.mean_gap_size) && !isnan(row.C_overlap_density), df_final)


# Add microscopic epsilon so Gamma distribution handles exact zeros
y_dens = Float64.(df_final.C_overlap_density) .+ 1e-6 
x_gaps = Float64.(df_final.mean_gap_size)
df_final.C_overlap_density_gamma = y_dens

# ---------------------------------------------------------
# STEP 1: The Sloping Single Hinge (The Physical Wall)
# ---------------------------------------------------------
println("\n[Step 1] Executing Sloping Hinge Optimization...")
best_sse = Inf
best_T_hinge = 0.0
best_hinge_mod = nothing

for T in 100.0:0.5:200.0
    df_final.regime_1 = df_final.mean_gap_size
    df_final.regime_2 = max.(0.0, T .- df_final.mean_gap_size)
    
    mod = lm(@formula(C_overlap_density ~ regime_1 + regime_2), df_final)
    preds = predict(mod)
    sse = sum((df_final.C_overlap_density .- preds).^2)
    
    if sse < best_sse
        best_sse = sse
        best_T_hinge = T
        best_hinge_mod = mod
    end
end

println("--> Primary Physical Wall (Hinge): $(best_T_hinge) bp")

# ---------------------------------------------------------
# STEP 2: 3-State Gamma Log-Likelihood (The Population Shifts)
# ---------------------------------------------------------
println("\n[Step 2] Executing 3-State Gamma Log-Likelihood Search...")
best_ll_3 = -Inf
best_T1 = 0.0
best_T2 = 0.0

for T1 in 120.0:1.0:145.0 
    for T2 in 80.0:1.0:(T1 - 5.0)
        y_relax = y_dens[x_gaps .>= T1]
        y_trans = y_dens[T2 .<= x_gaps .< T1]
        y_extrm = y_dens[x_gaps .< T2]
        
        if length(y_relax) > 20 && length(y_trans) > 20 && length(y_extrm) > 20
            ll = sum(logpdf.(fit(Gamma, y_relax), y_relax)) + 
                 sum(logpdf.(fit(Gamma, y_trans), y_trans)) + 
                 sum(logpdf.(fit(Gamma, y_extrm), y_extrm))
                 
            if ll > best_ll_3
                best_ll_3 = ll
                best_T1 = T1
                best_T2 = T2
            end
        end
    end
end

println("--> Population Transitional Shift (T1): $(best_T1) bp")
println("--> Population Extreme Rupture (T2):  $(best_T2) bp")

# ---------------------------------------------------------
# STEP 3: GMM Autopsy on the True Extreme Zone (< T2)
# ---------------------------------------------------------
println("\n[Step 3] Executing GMM Autopsy on Extreme Zone (< $(best_T2) bp)...")

df_extreme = filter(r -> r.mean_gap_size < best_T2, df_final)
extreme_vals = Float64.(df_extreme.C_overlap_density)

if nrow(df_extreme) < 10
    println("Not enough data in the Extreme zone for GMM.")
else
    # Reshape the 1D vector into an N x 1 matrix for the GMM package
    extreme_mat = reshape(extreme_vals, length(extreme_vals), 1)
    
    # Fit a 2-component Gaussian Mixture Model
    gmm = GMM(2, extreme_mat)
    
    # Extract probabilities
    ll, post = gmmposterior(gmm, extreme_mat)
    
    # Identify which column belongs to the "Higher Density" peak
    high_density_idx = gmm.μ[1] > gmm.μ[2] ? 1 : 2
    
    # Assign genomes based on >50% probability of belonging to the high-density peak
    probs_high = post[:, high_density_idx]
    df_extreme.extreme_type = [p > 0.5 ? "Hyper-Entangled" : "Efficient-Reduced" for p in probs_high]
    
    println("\n--> Biological Divergence in the Extreme Zone:")
    stats = combine(groupby(df_extreme, :extreme_type), 
        nrow => :count,
        :C_overlap_density => mean => :avg_overlap_density,
        :genome_size => mean => :avg_genome_size,
        :total_genes => mean => :avg_total_genes
    )
    println(stats)
end

println("\n======================================================")
println("Pipeline Complete. Ready for manuscript integration.")
println("======================================================")