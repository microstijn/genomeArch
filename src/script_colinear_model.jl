using DataFrames
using GLM
using Statistics
using Random
using LinearAlgebra

# ---------------------------------------------------------
# 1. VARIANCE INFLATION FACTOR (VIF)
# Diagnoses exactly which variables are dangerously covariant
# (VIF > 5 is concerning, VIF > 10 is severe collinearity)
# ---------------------------------------------------------
function calculate_vif(df::DataFrame, formula::FormulaTerm)
    # Extract the design matrix
    mf = ModelFrame(formula, df)
    mm = ModelMatrix(mf)
    X = mm.m
    names_X = coefnames(mf)
    
    n_vars = size(X, 2)
    vifs = Float64[]
    
    # Skip intercept (usually column 1)
    for i in 2:n_vars
        y_col = X[:, i]
        # X matrix without the i-th column
        X_other = X[:, setdiff(1:n_vars, i)]
        
        # Fit linear model: X_i ~ X_other
        beta = X_other \ y_col
        y_pred = X_other * beta
        
        # Calculate R^2
        ss_tot = sum((y_col .- mean(y_col)).^2)
        ss_res = sum((y_col .- y_pred).^2)
        r2_val = 1.0 - (ss_res / ss_tot)
        
        # VIF formula
        push!(vifs, 1.0 / (1.0 - r2_val))
    end
    
    return DataFrame(Predictor = names_X[2:end], VIF = vifs)
end

# ---------------------------------------------------------
# 2. K-FOLD CROSS-VALIDATION
# Proves the model isn't just memorizing the dataset
# ---------------------------------------------------------
function k_fold_cv(df::DataFrame, formula::FormulaTerm, target_col::Symbol; k=5, weights_col=nothing)
    n = nrow(df)
    indices = shuffle(1:n)
    fold_size = div(n, k)
    
    rmses = Float64[]
    
    for i in 1:k
        test_idx = indices[((i-1)*fold_size + 1):min(i*fold_size, n)]
        train_idx = setdiff(indices, test_idx)
        
        train_data = df[train_idx, :]
        test_data = df[test_idx, :]
        
        # Fit on K-1
        if isnothing(weights_col)
            mod = lm(formula, train_data)
        else
            mod = lm(formula, train_data, wts=train_data[!, weights_col])
        end
        
        # Predict on 1
        preds = predict(mod, test_data)
        actuals = test_data[!, target_col]
        
        # Calculate Root Mean Square Error
        rmse = sqrt(mean((preds .- actuals).^2))
        push!(rmses, rmse)
    end
    
    println("K-Fold RMSEs: ", round.(rmses, digits=2))
    println("Mean CV RMSE: ", round(mean(rmses), digits=2), " ± ", round(std(rmses), digits=2))
    return mean(rmses)
end

# ---------------------------------------------------------
# 3. BOOTSTRAPPING
# Bypasses parametric assumptions to give true 95% CIs
# ---------------------------------------------------------
function bootstrap_model(df::DataFrame, formula::FormulaTerm; iterations=1000, weights_col=nothing)
    n = nrow(df)
    
    # Fit base model to get coefficient names and count
    base_mf = ModelFrame(formula, df)
    coef_names = coefnames(base_mf)
    n_coef = length(coef_names)
    
    boot_coefs = zeros(iterations, n_coef)
    
    for i in 1:iterations
        # Sample with replacement
        boot_idx = rand(1:n, n)
        boot_df = df[boot_idx, :]
        
        if isnothing(weights_col)
            mod = lm(formula, boot_df)
        else
            mod = lm(formula, boot_df, wts=boot_df[!, weights_col])
        end
        
        boot_coefs[i, :] = coef(mod)
    end
    
    # Calculate 2.5% and 97.5% quantiles for 95% CI
    results = DataFrame(
        Predictor = coef_names,
        Mean_Estimate = [mean(boot_coefs[:, j]) for j in 1:n_coef],
        Lower_95 = [quantile(boot_coefs[:, j], 0.025) for j in 1:n_coef],
        Upper_95 = [quantile(boot_coefs[:, j], 0.975) for j in 1:n_coef]
    )
    
    return results
end

# ---------------------------------------------------------
# EXECUTION
# ---------------------------------------------------------
# Assuming 'df_clean' is the output from your optimize_models function
formula_6 = @formula(C_overlap_nr ~ genome_size + OGT + genome_size & lifestyle_tier + compression_penalty + compression_penalty_sq + compression_penalty_sq & lifestyle_tier)
formula_6 = @formula(C_overlap_nr ~ genome_size + mean_gap_size + compression_penalty + compression_penalty_sq + genome_size & lifestyle_tier)

formula_ultimate = @formula(C_overlap_nr ~ genome_size + compression_penalty + compression_penalty_sq + genome_size & lifestyle_tier)

df_temp = processed_df

df_temp.overlap_prob = df_temp.C_overlap_nr ./ df_temp.total_genes
df_temp.compression_penalty = max.(0.0, threshold .- df_temp.mean_gap_size)
df_temp.compression_penalty_sq = max.(0.0, threshold .- df_temp.mean_gap_size) .^ 2
    
println("\n--- 1. Collinearity Check ---")
calculate_vif(df_temp, formula_6)

println("\n--- 2. Out-of-Sample Predictive Power ---")
k_fold_cv(df_temp, formula_6, :C_overlap_nr, k=10, weights_col=:model_weights)

println("\n--- 3. Robust Bootstrapped Confidence Intervals ---")
display(bootstrap_model(df_temp, formula_6, iterations=1000, weights_col=:model_weights))