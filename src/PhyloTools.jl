module PhyloTools

using Phylo
using PhyloNetworks # Added for robust Newick parsing
using DataFrames
using CSV

export prune_gtdb_tree, inspect_tree_file

"""
    inspect_tree_file(gtdb_tree_file::String)

Reads and parses a Newick tree file using the robust PhyloNetworks parser, 
then prints a summary of its structure.
"""
function inspect_tree_file(gtdb_tree_file::String)
    println("Inspecting tree file: ", gtdb_tree_file)
    
    if !isfile(gtdb_tree_file)
        @error "Tree file not found at '$gtdb_tree_file'."
        return
    end
    
    # --- 1. Parse the Tree using PhyloNetworks.jl ---
    println("Parsing the tree with PhyloNetworks.jl...")
    # readTopology is more robust and handles complex Newick formats like GTDB's.
    tree = readTopology(gtdb_tree_file)
    println("Tree parsed successfully.")
    
    # --- 2. Print a Summary of the Tree ---
    # Note: PhyloNetworks objects have a different API than Phylo.jl objects
    println("\n--- Tree Summary (from PhyloNetworks) ---")
    println("Number of tips (leaves): ", length(tree.leaf))
    println("Number of internal nodes: ", length(tree.node) - length(tree.leaf))
    println("Total nodes: ", length(tree.node))
    
    println("\nExample tip names (first 5):")
    for i in 1:min(5, length(tree.leaf))
        println("  - ", tree.leaf[i].name)
    end
    println("----------------------------------------")

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

    # --- 2. Load the full GTDB tree using the robust PhyloNetworks parser ---
    println("Loading the full GTDB tree from: ", gtdb_tree_file)
    if !isfile(gtdb_tree_file)
        @error """
        GTDB tree file not found at '$gtdb_tree_file'.
        Please download the latest tree from the GTDB website (e.g., bac120_sp_r220.tree or ar53_r220.tree)
        """
        return
    end
    
    full_tree = readTopology(gtdb_tree_file)
    println("Tree parsed successfully with PhyloNetworks.")

    # --- 3. Prune the tree ---
    println("Pruning the tree to keep only the specified accessions...")
    # The function `deleteleaf!` is the standard way to prune in PhyloNetworks.
    tips_to_delete = setdiff(Set(n.name for n in full_tree.leaf), accessions_to_keep)
    
    if length(tips_to_delete) == length(full_tree.leaf)
        @error "None of the provided accessions were found in the GTDB tree. Check accession format (e.g., 'RS_GCF_...')."
        return
    end
    
    for tip in tips_to_delete
        deleteleaf!(full_tree, tip)
    end

    # --- 4. Write the pruned tree to a new file ---
    println("Writing pruned tree to: ", output_tree_file)
    mkpath(dirname(output_tree_file))
    # `writeTopology` is the function for writing PhyloNetworks trees.
    writeTopology(full_tree, output_tree_file)

    println("\nTree pruning complete. Pruned tree contains $(length(full_tree.leaf)) tips.")
end

end # end module

