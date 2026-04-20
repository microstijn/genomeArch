using DataFrames
using CSV

# ==============================================================================
# 1. LOAD DATASETS
# ==============================================================================
println("Loading datasets for dynamic merging...")

path_master = raw"D:\pipeline_output\merged_imputed_ogt.csv"
path_madin  = raw"D:\bacteria_archaea_traits.csv"  # <-- Update this path

df_master = CSV.read(path_master, DataFrame)
df_madin  = CSV.read(path_madin, DataFrame, missingstring = "NA")


original_row_count = nrow(df_master)
println("Target Row Count to Maintain: $original_row_count")

df_master.taxId = string.(df_master.taxId)
df_madin.tax_id = string.(df_madin.tax_id)
df_madin.species_tax_id = string.(df_madin.species_tax_id)
# Filter overlapping columns
overlapping_cols = intersect(names(df_master), names(df_madin))
cols_to_keep = setdiff(names(df_madin), overlapping_cols)
push!(cols_to_keep, "tax_id")
madin_subset = select(df_madin, unique(cols_to_keep))

println("Dereplicating Madin to mathematically prevent row explosion...")
gdf = groupby(madin_subset, :tax_id)
agg_cols = setdiff(names(madin_subset), ["tax_id"])

madin_derep = combine(gdf, agg_cols .=> (col -> begin
    valid_vals = collect(skipmissing(col))
    if isempty(valid_vals) missing
    elseif eltype(valid_vals) <: Number median(valid_vals)
    else first(valid_vals) end
end) .=> agg_cols)

# ==============================================================================
# 3. MERGE & VERIFY
# ==============================================================================
df_merged = leftjoin(df_master, madin_derep, on = :taxId => :tax_id)

if nrow(df_merged) != original_row_count
    error("CRITICAL MERGE FAILURE: Rows changed from $original_row_count to $(nrow(df_merged))")
else
    println("SUCCESS: Merge completed. Row count perfectly matches original ($original_row_count rows).")
end

# ==============================================================================
# 4. RUN TAXONOMIC IMPUTATION
# ==============================================================================
println("Loading Taxonomy DB...")
# Ensure you point this to your actual NCBI Taxonomy dump locations

nodesTax = raw"D:\ncbi_downloads\taxdump\nodes.dmp"
namesTax = raw"D:\ncbi_downloads\taxdump\names.dmp"
db = Taxonomy.DB(nodesTax, namesTax) 

# Define exactly what we want to impute, and the mathematical rule to use
traits_to_impute = [
    (:doubling_h, median),       # Continuous -> Median
    (:gc_content, median),       # Continuous -> Median
    (:rRNA16S_genes, median),    # Continuous -> Median
    (:gram_stain, get_mode),     # Categorical -> Mode
    (:motility, get_mode)        # Categorical -> Mode
]

for (trait, func) in traits_to_impute
    if string(trait) in names(df_merged)
        println("--> Imputing $trait...")
        # Re-assign df_merged to accumulate the new imputed columns safely
        global df_merged = impute_by_taxonomy(df_merged, trait, :taxId, db, agg_func=func)
    else
        println("Warning: Column $trait not found, skipping imputation.")
    end
end

# ==============================================================================
# 5. SAVE FINAL DATASET
# ==============================================================================
output_path = raw"D:\pipeline_output\master_imputed_genarch_with_madin.csv"
CSV.write(output_path, df_merged)
println("Final dataset saved to: ", output_path)
println("Final Row Count: $(nrow(df_merged))")

sum(ismissing.(df_merged.doubling_h))
