module TaxonomyTools

export process_taxonomy

using CSV
using DataFrames
using JSON3
using Taxonomy

"""
    get_nested(obj, keys...)

Safely access nested keys in a JSON3 object. Returns `missing` if any key is not found.
"""
function get_nested(obj, keys...)
    val = obj
    for key in keys
        if !haskey(val, key)
            return missing
        end
        val = val[key]
    end
    return val
end

"""
    getFullTaxonomy(taxId, taxNodesLoc, taxNamesLoc)

Robustly retrieves the full taxonomic lineage for a vector of NCBI TaxIDs.
Handles missing ranks by filling them with "NA".
"""
function getFullTaxonomy(taxId, taxNodesLoc::String, taxNamesLoc::String)
    
    errNode = isfile(taxNodesLoc)
    errNode == true || throw("$taxNodesLoc is not a file")
    errName = isfile(taxNamesLoc)
    errName == true || throw("$taxNamesLoc is not a file")

    println("Loading taxonomy database from taxdump files...")
    db = Taxonomy.DB(taxNodesLoc, taxNamesLoc)

    taxid_ordered_dataframe = DataFrame(
        superkingdom = String[],
        phylum = String[],
        class = String[],
        order = String[],
        family = String[],
        genus = String[],
        species = String[]
    )

    total_ids = length(taxId)
    println("Looking up taxonomy for $total_ids entries...")
    for (i, id) in enumerate(taxId)
        # Set defaults
        superkingdom, phylum, class, order, family, genus_tax, species_tax = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
        
        if ismissing(id)
            push!(taxid_ordered_dataframe, ["NA", "NA", "NA", "NA", "NA", "NA", "NA"])
            continue
        end

        lineage_tax = "NA"
        try 
            tax = Taxon(id, db)
            lineage_tax = Lineage(tax)
        catch;
            lineage_tax = "NA"
        end

        if lineage_tax != "NA"
            try 
                superkingdom = string(lineage_tax[:superkingdom])
            catch; # Already "NA"
            end

            try 
                phylum = string(lineage_tax[:phylum])
            catch; # Already "NA"
            end

            try 
                class = string(lineage_tax[:class])
            catch; # Already "NA"
            end

            try 
                order = string(lineage_tax[:order])
            catch; # Already "NA"
            end

            try 
                family = string(lineage_tax[:family])
            catch; # Already "NA"
            end
    
            try 
                genus_tax = string(lineage_tax[:genus])
            catch; # Already "NA"
            end

            try 
                species_tax = string(lineage_tax[:species])
            catch; # Already "NA"
            end
        end

        push!(taxid_ordered_dataframe, [
            superkingdom,
            phylum,
            class,
            order,
            family,
            genus_tax,
            species_tax
        ])
    end
    println("\nTaxonomy lookup complete.")
    return taxid_ordered_dataframe
end


"""
    process_taxonomy(jsonl_files, taxdump_dir, output_dir)

Main exported function. Reads NCBI JSONL files, extracts metadata,
retrieves full taxonomy, and writes the combined data to a CSV file.
"""
function process_taxonomy(jsonl_files::Vector{String}, taxdump_dir::String, output_dir::String)
    nodes_dmp = joinpath(taxdump_dir, "nodes.dmp")
    names_dmp = joinpath(taxdump_dir, "names.dmp")

    for jsonl_path in jsonl_files
        if !isfile(jsonl_path)
            @warn "File not found: $jsonl_path. Skipping."
            continue
        end

        println("Processing file: $jsonl_path")
        
        # Pre-allocate vectors
        organismName = Union{Missing, String}[]
        taxId = Union{Missing, Int64}[]
        checkmSpeciesTaxId = Union{Missing, Int}[]
        accession = Union{Missing, String}[]
        completeness = Union{Missing, Float64}[]
        gcPercent = Union{Missing, Float64}[]
        
        open(jsonl_path, "r") do io
            for line in eachline(io)
                record = JSON3.read(line)
                
                push!(organismName, get_nested(record, :organism, :organismName))
                push!(taxId, get_nested(record, :organism, :taxId))
                push!(checkmSpeciesTaxId, get_nested(record, :checkmInfo, :checkmSpeciesTaxId))
                push!(accession, get_nested(record, :accession))
                push!(completeness, get_nested(record, :checkmInfo, :completeness))
                push!(gcPercent, get_nested(record, :assemblyStats, :gcPercent))
            end
        end
        println("JSONL processing complete.")

        # Create DataFrame from collected data
        df = DataFrame(
            organismName = organismName,
            taxId = taxId,
            checkmSpeciesTaxId = checkmSpeciesTaxId,
            accession = accession,
            completeness = completeness,
            gcPercent = gcPercent
        )

        full_tax = getFullTaxonomy(df.taxId, nodes_dmp, names_dmp)
        
        # Combine the dataframes
        df_final = hcat(df, full_tax)

        # Write to output file
        output_filename = first(split(basename(jsonl_path), '.')) * "_TaxId.csv"
        output_path = normpath(joinpath(output_dir, output_filename))
        
        println("Writing results to: $output_path")
        try
            CSV.write(output_path, df_final)
        catch e
            @error "Failed to write CSV file at $output_path." exception=(e, catch_backtrace())
        end
        println("Finished processing $jsonl_path.\n")
    end
end

end # end of module