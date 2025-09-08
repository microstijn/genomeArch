module genArch

# Include the code for each module
include("TaxonomyTools.jl")
include("EnvironmentTools.jl")
include("ArchitectureTools.jl")
include("TestFunctions.jl")
include("AnalysisTools.jl")
include("PhyloTools.jl")


# Make the main functions from each module available to the user
using .TaxonomyTools
export process_taxonomy

using .EnvironmentTools
export fetch_environments

using .ArchitectureTools
export calculate_architecture

using .TestFunctions
export run_all_tests

using .AnalysisTools
export consolidate_to_genomes
export calculate_density_score
export perform_pic_analysis

using .PhyloTools
export prune_gtdb_tree
export inspect_tree_file

end # module