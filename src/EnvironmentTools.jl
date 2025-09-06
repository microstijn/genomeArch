module EnvironmentTools

export fetch_environments

using CSV
using DataFrames
using HTTP
using JSON3
using ProgressMeter

#-----------------------------------------------------------------
# Helper Function for API call
#-----------------------------------------------------------------
function fetch_omnicrobe_env(taxid::Int, base_url::String)
    full_url = string(base_url, taxid)
    try
        r = HTTP.get(full_url, status_exception=false)

        if r.status == 200
            json_obj = JSON3.read(r.body)
            if isempty(json_obj)
                return (taxid, missing, missing)
            end
            
            envs = [item.obt_forms[1] for item in json_obj]
            obtids = [item.obtid for item in json_obj]
            return (taxid, envs, obtids)
        else
            if r.status != 404
                @warn "Request failed for taxid $taxid with status $(r.status)"
            end
            return (taxid, missing, missing)
        end
    catch e
        @error "Exception for taxid $taxid: $e"
        return (taxid, missing, missing)
    end
end

#-----------------------------------------------------------------
# Main Exported Function
#-----------------------------------------------------------------
function fetch_environments(input_csv::String, output_dir::String)
    if !isfile(input_csv)
        @error "Input file not found: $input_csv"
        return
    end

    mkpath(output_dir)
    
    # --- FILENAME CORRECTION IS HERE ---
    # The output files will now correctly have a .tsv extension.
    base_name = first(split(basename(input_csv), '_'))
    genome_report_path = joinpath(output_dir, base_name * "_genomes_report.tsv")
    environments_path = joinpath(output_dir, base_name * "_genome_environments.tsv")

    println("Reading input file: $input_csv")
    df = CSV.File(input_csv) |> DataFrame

    base_url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A"

    unique_taxids = unique(skipmissing(df.taxId))
    n_taxids = length(unique_taxids)

    println("Querying OmniMicrobe for $n_taxids unique taxon IDs...")
    
    # Use asyncmap for a clean, concurrent implementation with a progress bar
    results = @showprogress "Querying API..." asyncmap(taxid -> fetch_omnicrobe_env(taxid, base_url), unique_taxids)

    println("\nAPI queries complete. Processing results...")

    taxid_to_accession = Dict(zip(df.taxId, df.accession))

    env_df = DataFrame(accession=String[], environment=String[], obtId=String[])
    
    for (taxid, envs, obtids) in results
        if !ismissing(envs)
            accession = get(taxid_to_accession, taxid, "Unknown_Accession")
            for (env, obtid) in zip(envs, obtids)
                push!(env_df, (accession, env, obtid))
            end
        end
    end

    println("Writing tidy environment data to: $environments_path")
    CSV.write(environments_path, env_df, delim='\t')
    
    if "environments" in names(df)
        select!(df, Not(:environments))
    end
    if "obtId" in names(df)
        select!(df, Not(:obtId))
    end
    
    println("Writing main genome report to: $genome_report_path")
    CSV.write(genome_report_path, df, delim='\t')

    println("Environment fetching complete.")
end

end # end of module