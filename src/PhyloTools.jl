module PhyloTools

using Phylo
using DataFrames
using CSV

export prune_gtdb_tree, inspect_tree_file

"""
    inspect_tree_file(gtdb_tree_file::String)

Reads, cleans, and parses a Newick tree file, then prints a summary of its structure.
"""
function inspect_tree_file(gtdb_tree_file::String)
    println("Inspecting tree file: ", gtdb_tree_file)
    
    if !isfile(gtdb_tree_file)
        @error "Tree file not found at '$gtdb_tree_file'."
        return
    end
    
    # --- 1. Read and Clean the Tree String ---
    println("Reading and cleaning tree string...")
    tree_string = read(gtdb_tree_file, String)
    # This regex replicates the functionality of the GTDB-Tk `tree_to_itol.py` script,
    # removing the internal node labels (bootstraps) that cause parsing errors.
    cleaned_tree_string = replace(tree_string, r"\)([^,:]*):" => "):")
    
    # --- 2. Parse the Cleaned String ---
    println("Parsing the cleaned tree...")
    tree = parsenewick(cleaned_tree_string)
    println("Tree parsed successfully.")
    
    # --- 3. Print a Summary of the Tree ---
    leaf_names = getleafnames(tree)
    num_leaves = length(leaf_names)
    num_nodes = nnodes(tree)

    println("\n--- Tree Summary ---")
    println("Total number of nodes: ", num_nodes)
    println("Number of tips (leaves): ", num_leaves)
    println("Number of internal nodes: ", num_nodes - num_leaves)
    
    println("\nExample tip names (first 5):")
    for i in 1:min(5, num_leaves)
        println("  - ", leaf_names[i])
    end
    println("--------------------")

    return tree
end


"""
    prune_gtdb_tree(gtdb_tree_file::String, accessions_file::String, output_tree_file::String)

Prunes a full GTDB reference tree to only include the accessions present in the user's dataset.
"""
function prune_gtdb_tree(gtdb_tree_file::String, accessions_file::String, output_tree_file::String)
    println("Starting GTDB tree pruning...")
    
    # --- 1. Load the list of accessions to keep ---
    println("Reading accessions from: ", accessions_file)
    df = CSV.File(accessions_file) |> DataFrame
    # GTDB accessions are in the format 'RS_GCF_...' or 'GB_GCA_...', so we need to reformat ours.
    accessions_to_keep = Set("RS_" .* df.genome_name)
    union!(accessions_to_keep, Set("GB_" .* df.genome_name))
    println("Found $(length(df.genome_name)) unique genome accessions to keep.")

    # --- 2. Load and Clean the full GTDB tree string ---
    println("Loading and cleaning the full GTDB tree from: ", gtdb_tree_file)
    if !isfile(gtdb_tree_file)
        @error """
        GTDB tree file not found at '$gtdb_tree_file'.
        Please download the latest tree from the GTDB website (e.g., bac120_sp_r220.tree or ar53_r220.tree)
        """
        return
    end
    
    tree_string = read(gtdb_tree_file, String)
    cleaned_tree_string = replace(tree_string, r"\)([^,:]*):" => "):")
    
    full_tree = parsenewick(cleaned_tree_string)

    # --- 3. Prune the tree ---
    println("Pruning the tree to keep only the specified accessions...")
    tips_in_tree = Set(getleafnames(full_tree))
    final_tips = intersect(accessions_to_keep, tips_in_tree)
    
    if isempty(final_tips)
        @error "None of the provided accessions were found in the GTDB tree. Check accession format (e.g., 'RS_GCF_...' or 'GB_GCA_...')."
        return
    end
    
    keeponly!(full_tree, final_tips)

    # --- 4. Write the pruned tree to a new file ---
    println("Writing pruned tree to: ", output_tree_file)
    mkpath(dirname(output_tree_file))
    open(output_tree_file, "w") do io
        print(io, full_tree)
    end

    println("\nTree pruning complete. Pruned tree contains $(length(final_tips)) tips.")
end

end # end module

