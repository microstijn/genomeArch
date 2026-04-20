using CSV
using DataFrames

# Load your data
df = CSV.read(raw"D:\pipeline_output\overlapping_pairs_annotated.csv", DataFrame)

outfold = raw"D:\pipeline_output"

# Isolate the survivors
opp_strand_df = filter(row -> row.overlap_type in ["Convergent", "Divergent"], df)
survivors_df = filter(row -> row.overlap_length > 120, opp_strand_df)

# Get the unique genomes
unique_genomes = unique(survivors_df.genome)

# Write to a text file for the NCBI downloader
open(joinpath(outfold, "target_genomes.txt"), "w") do io
    for g in unique_genomes
        println(io, g)
    end
end

# Save the subset so our next script can loop through it fast
CSV.write(joinpath(outfold, "survivor_overlaps.csv"), survivors_df)
println("Saved $(length(unique_genomes)) genomes to target_genomes.txt")