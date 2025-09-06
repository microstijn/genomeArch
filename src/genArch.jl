module genArch

# Include the code for each module
include("TaxonomyTools.jl")
include("EnvironmentTools.jl")
include("ArchitectureTools.jl")

# Make the main functions from each module available to the user
using .TaxonomyTools
export process_taxonomy

using .EnvironmentTools
export fetch_environments

using .ArchitectureTools
export calculate_architecture

end # module