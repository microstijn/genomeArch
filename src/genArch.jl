module genArch

# Include the code for each module
include("TaxonomyTools.jl")
include("EnvironmentTools.jl")
include("ArchitectureTools.jl")
include("TestFunctions.jl")
include("AnalysisTools.jl")


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

end # module