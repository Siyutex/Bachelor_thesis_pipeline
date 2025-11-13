# generate a phylogenetic tree from X_CNV (so the cnvs are the evolutionary steps)
# show the tree
# save in format where we can select clades corresponding to cells in adata so we can define transition cells

import scanpy as sc
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj
import helper_functions as hf
import os
from Bio import Phylo
import warnings
import copy

def build_tree(adata, cnv_score_matrix, distance_metric):# get cnv score matrix

    vprint("getting cnv score matrix")
    X = adata.obsm[cnv_score_matrix].toarray() # need to turn scipy sparse csr matrix into numpy array for pdist to work

    # DEBUGGING, print percentage of NaNs in X
    vprint(f"Percentage of NaNs in cnv score matrix: {np.isnan(X).sum() / X.size * 100}%")
    vprint(f"Number of NaNs in cnv score matrix: {np.isnan(X).sum()}")

    # save cell identifiers to list
    labels = adata.obs_names.tolist()

    # add a tiny jitter so pearson correlation can be computed for all cells 
    # (if there are rows with 0 variance (ie all 0s) then pearson correlation cannot be computed --> NaNs)
    epsilon = 1e-16
    vprint("adding jitter to X")
    X_jittered = X + epsilon * np.random.randn(*X.shape)

    # Compute distance matrix
    print("computing distance matrix")
    D = squareform(pdist(X_jittered, metric=distance_metric)) # correlation = how similar cnv profiles are, invariant to magnitude of expression
    vprint(f"Percentage of NaNs in distance matrix: {np.isnan(D).sum() / D.size * 100}%")
    vprint(f"Number of NaNs in distance matrix: {np.isnan(D).sum()}")
    dm = DistanceMatrix(D, ids=labels)

    # Build neighbor-joining tree
    print("building neighbor-joining tree")
    tree = nj(dm) # this algorithm is deterministic
    return tree

def divide_tree(tree, n_clades: int) -> list[tuple[str, ...]]: 

    def all_parents(tree):
        parents = {}
        for clade in tree.find_clades(order="level"):
            for child in clade:
                parents[child.name] = clade
        return parents # dict key = node name, value = parent (clade object)

    # visitation order (widest clades first)
    def get_visitation_order_width(tree):
        visitation_order = {}
        for clade in tree.find_clades(order="level"):
            children = list(clade.clades)
            non_terminal_children = list(set(children).difference(set(terminals)))

            width_list = []
            for ntc in non_terminal_children:
                width = len(list(ntc.get_terminals())) # width of the clade
                width_list.append((ntc, width)) # tuple of child as clade object and its width
            
            width_list.sort(key=lambda x: x[1], reverse=True) # sort in descending order by width
            visitation_order[clade.name] = width_list
        
        return visitation_order
    
    
    # return the preterminal node reached if you always take the first child in the passed visitation order
    def find_start_node(tree, visitation_order):
        start_node = None
        
        def recurse(node):
            nonlocal start_node
            if all(child in terminals for child in node.clades):
                start_node = node # set start node to the first node encountered that only has terminal children (ie the preterminal node)
                return

            next_node = visitation_order[node.name][0][0] # first tuple element in first tuple = clade object

            recurse(next_node)

        recurse(tree.root)

        return start_node
    

    # recursive function to build clades
    def recurse(node, child_n_leaves):
        leaves = node.get_terminals() # leaves of current clade
        n_leaves = len(leaves) # size of current clade
        nonlocal clade_root


        if any(leaf not in leaf_pool for leaf in leaves): # if any of the current node's leaves already belong to another clade
            return False # this node cannot be chosen, its child must be chosen instead (if we choose this node, then overlapping clades happen)

        if n_leaves >= (target_clade_size):
            if abs(child_n_leaves - target_clade_size) < abs(n_leaves - target_clade_size): # if child is closer choose child
                return False # false means this node is not chosen
            else:
                clade_leaf_sets.append(tuple(leaves)) # add leaves to clade (clade objects)
                for leaf in clade_leaf_sets[-1]:
                    leaf_pool.remove(leaf) # remove matching clade objects from leaf pool
                clade_leaf_sets[-1] = tuple(leaf.name for leaf in clade_leaf_sets[-1]) # change clade objects to strings of the leaf name for downstream processing
                clade_root = node # set nonlocal variable "clade_root" to chosen node
                return True # true means this node is chosen
        
        if parents[node.name] != internal_tree.root:
            parent_chosen = recurse(parents[node.name], n_leaves) # recursive call on parent
        else:
            print("recursion reached root node of tree, interrupting recursion")
            print(f"the tree root's child {node.name} was chosen as the calde root instead")
            parent_chosen = False # the root cannot be chosen since that would include the entire tree so its child must be chosen instead

        if parent_chosen == False: # if the parent was not chosen (because this node is closer to ideal clade size)
            clade_leaf_sets.append(tuple(leaves)) # add leaves to clade (clade objects)
            for leaf in clade_leaf_sets[-1]:
                leaf_pool.remove(leaf) # remove matching clade objects from leaf pool
            clade_leaf_sets[-1] = tuple(leaf.name for leaf in clade_leaf_sets[-1]) # change clade objects to strings of the leaf name for downstream processing
            clade_root = node # set nonlocal clade root to chosen node
        return None # None means this node OR one of its ancestors was chosen
    

    # sanity check
    if n_clades < 1:
        raise ValueError("n_clades must be at least 1")
    if n_clades == 1:
        return [tuple(tip.name for tip in internal_tree.tips())]
    
    # define internal tree to not change original
    internal_tree = copy.deepcopy(tree)
    
    # local vars
    terminals = list(internal_tree.get_terminals())
    parents = all_parents(internal_tree)
    leaf_pool = list(internal_tree.get_terminals(order = "postorder")) # list of available leaves (names of the cells) for selection
    n_leaves = len(leaf_pool) # total number of leaves in tree
    target_clade_size = n_leaves / n_clades # target number of leaves in ideal clade
    clade_leaf_sets = [] # list of chosen clades (list of str tuples)

    # create clades
    i = 0
    while len(clade_leaf_sets) < n_clades and len(leaf_pool) > 1: # repeat until n_clades reached or no leaves left
        print(f"Building clade {str(i)}")
        visitation_order = get_visitation_order_width(internal_tree) # compute visitation order for each run (so nodes in existing clades cannot be visited)
        start_node = find_start_node(internal_tree, visitation_order) # find start node
        clade_root = None # root of clade (will be set in recurse)
        recurse(start_node, 0) # recurse upward from that node and add entry to clade_leaf_sets
        vprint(f"Clade root for clade {i}: {clade_root.name}") # debugging info
        parents[clade_root.name].clades.remove(clade_root) # remove clade root and its descendants from tree
        i += 1

        if len(clade_leaf_sets) >= n_clades:
            vprint("Interrupting loop, reached target number of clades.")
            vprint("Number of leaves remaining: " + str(len(leaf_pool)))
        elif len(leaf_pool) <= 1:
            vprint("Interrupting loop, all leaves used.")
            vprint("Number of clades built: " + str(len(clade_leaf_sets)))



    # info about clades
    vprint(f"Number of clades: {len(clade_leaf_sets)}")
    for i, leaves in enumerate(clade_leaf_sets):
        vprint(f"Clade {i}: {len(list(leaves))}")

    return clade_leaf_sets


def fix_clades(tree, clade_list):
    # clade list is a list of string tuples (one tuple per clade, one string per cell ID in the clade)

    # same visitation order as in divide tree (where clades are built)
    def get_visitation_order_width(tree):
        # key = node.name value = list of tuples of child as clade object and its width
        # ever node in the tree has a key in this dict, even terminals (but their value is [])
        visitation_order = {}
        for clade in tree.find_clades(order="level"):
            children = list(clade.clades)

            width_list = []
            for child in children:
                width = len(list(child.get_terminals())) # width of the clade
                width_list.append((child, width)) # tuple of child as clade object and its width
            
            width_list.sort(key=lambda x: x[1], reverse=True) # sort in descending order by width
            visitation_order[clade.name] = width_list
        
        return visitation_order
    

    def get_leaf_sequence(tree, visitation_order):
        # return a list with strings of leaf names ordered by visitation order (so leaf that is visited first is 1, then 2 ...)
        leaf_sequence = []

        def recurse(node):
            if node.clades == []: # append terminal nodes to the list
                leaf_sequence.append(node.name)
                return

            for i,_ in enumerate(node.clades): # recurse 
                next_node = visitation_order[node.name][i][0]
                recurse(next_node)

        recurse(tree.root)
        
        return leaf_sequence
    

    def get_fixed_clade_list(clade_list, leaf_sequence):

        def find_gaps(clade): # return a list of lists with cells at their proper indices according to leaf sequence, empty indices filled with 0
            # find cell with lowest index in leaf_sequence
            start_index = leaf_sequence.index(clade[0])
            for cell in clade:
                if cell in leaf_sequence and leaf_sequence.index(cell) < start_index:
                    start_index = leaf_sequence.index(cell)
            vprint(f"start index: {start_index}")
            # find cell with highest index in leaf_sequence
            end_index = leaf_sequence.index(clade[-1])
            for cell in clade:
                if cell in leaf_sequence and leaf_sequence.index(cell) > end_index:
                    end_index = leaf_sequence.index(cell)
            vprint(f"end index: {end_index}")
            # create a list of that length
            fixed_clade = [0 for _ in range(end_index - start_index + 1)] # (eg start = 0 end is 10, diff = 10, but 11 values)
            vprint(f"fixed clade length: {len(fixed_clade)} (including 0s for gaps)")
            vprint(f"clade length: {len(clade)}")
            # fill list with clade names at correct index
            for cell in clade:
                try:
                    index = leaf_sequence.index(cell) - start_index
                    fixed_clade[index] = cell
                except Exception as e:
                    print(f"Cell {cell} has index {leaf_sequence.index(cell)} in leaf sequence but fixed_clade has length {len(fixed_clade)}")
                    raise e
            return fixed_clade

        gapped_clade_list = []
        for clade in clade_list:
            gapped_clade = find_gaps(clade)
            gapped_clade_list.append(gapped_clade)

        # fixed_clade_list is now a list of lists, where each clade is a list of the string names of its leaf nodes, in the order at which they appear in the tree from back to front. If there is a gap in the sequence (ie another clade intersperses it), then that gap is filled with 0s)
        # so now we just find the 0 gaps in each fixed clade and split the clade at each gap


        fixed_clade_list = [] # this shall be a list where clades that have a gap are divided into n + 1 clades where n is the number of gaps
        for clade in gapped_clade_list:
            gap_start_indices = [] # first index of a gap
            for i,entry in enumerate(clade):
                if entry == 0 and clade[i-1] != 0  and i != 0: # make sure we are not checkoing the first index or it will wrap around to the last element
                    gap_start_indices.append(i)

            gap_end_indices = [] # last index of a gap
            for i,entry in enumerate(clade):
                if entry == 0 and clade[i+1] != 0 and i != len(clade)-1:
                    gap_end_indices.append(i)

            if gap_start_indices == [] and gap_end_indices == []: # if there are no gaps, no need to change the clade, just append as is
                vprint(f"clade {gapped_clade_list.index(clade)} has no gaps, appending as is")
                fixed_clade_list.append(clade)
            else:
                gaps = list(zip(gap_start_indices, gap_end_indices)) # list of tuples, each tuple contains the first and last index of a gap
                vprint(f"gaps in clade {gapped_clade_list.index(clade)} are: {gaps}")
                
                index_list = [] # index list should contain the indeces of where proper cell names start and end in pairs of 2 -> merge to tuple
                index_list.append(0)
                for gap in gaps:
                    index_list.append(gap[0]-1)
                    index_list.append(gap[1]+1)
                index_list.append(len(clade)-1) 
                
                filled_spaces = [[index_list[i], index_list[i+1]] for i in range(0, len(index_list)-1, 2)] # tuples of start and end indices of filled spaces (everything the is between the gaps, ie the cells)
                vprint(f"filled spaces in clade {gapped_clade_list.index(clade)} are: {filled_spaces}")
                for i,fs in enumerate(filled_spaces):
                    fixed_clade_list.append(clade[fs[0]:fs[1]])

        return fixed_clade_list
               
    internal_tree = copy.deepcopy(tree)
    visitation_order = get_visitation_order_width(internal_tree)
    leaf_sequence = get_leaf_sequence(internal_tree, visitation_order)
    fixed_clade_list = get_fixed_clade_list(clade_list, leaf_sequence)

    vprint(f"initial number of clades was: {len(clade_list)}")
    vprint(f"final number of clades (clades split at gaps) is: {len(fixed_clade_list)}")
    vprint(f"Clade lengths are {[len(clade) for clade in fixed_clade_list]}")

    return fixed_clade_list


def get_states(adata, metric, n_transition_clades, cutoff): 
    # metric is one of adata.obs["cancer_state"] or adata.obs["cancer_state_inferred"] (tells you which cells are cancerous and which non_cancerous)
    # cutoff is eg: calde has 90% c 10% nc -> classified as c, 89% c 11% nc, classified as unsure 
    


    vprint("Importing clade annotations")
    vprint(f"Number of NaNs in cnv clade: {np.sum(np.isnan(adata.obs['cnv_clade']))}")

    metrics = np.asarray(adata.obs[metric]) # array of metric column 
    clades = np.asarray(adata.obs["cnv_clade"])# array of clade column
    unique_clades = np.unique(clades) # array of unique clades
    vprint(f"Unique clades found: {unique_clades} ")

    entropy_dict = {}
    frac_c_dict = {}
    frac_nc_dict = {}

    # compute entropy for each clade
    for clade in unique_clades: 
        if clade != -1.: # -1 is not a clade it is just all cells that couldn't be assgined to a clade
            mask = clades == clade
            clade_metrics = metrics[mask]

            # get fractions of cancerous and non_cancerous cells in the clade
            n_total = len(clade_metrics)
            n_c = np.sum(clade_metrics == "cancerous")
            n_nc = np.sum(clade_metrics == "non_cancerous")

            if n_total == 0: # if no cells in clade
                frac_c = 0.0
                frac_nc = 0.0
                entropy = 0.0
                warnings.warn(f"No cells in clade {clade}, cannot compute entropy")
            else:
                frac_c = n_c / n_total
                frac_nc = n_nc / n_total
                # Shannon entropy
                if frac_c in (0.0, 1.0): # if clade is all cancerous or all non_cancerous then entropy is 0
                    entropy = 0.0
                else:
                    entropy = -(frac_c * np.log2(frac_c) + frac_nc * np.log2(frac_nc))

            # assign values in dictionary for current clade
            entropy_dict[clade] = entropy
            frac_c_dict[clade] = frac_c
            frac_nc_dict[clade] = frac_nc

    vprint("Clade summary:")
    for clade in unique_clades:
        if clade != -1.:
            vprint(f"Clade: {clade}, has {len(metrics[clades == clade])} cells. Cancerous fraction: {frac_c_dict[clade]}, Non-cancerous fraction: {frac_nc_dict[clade]}, Entropy: {entropy_dict[clade]}")

    # choose clades with highest entropy
    sorted_clades = sorted(entropy_dict, key=entropy_dict.get, reverse=True)
    transition_clades = set(sorted_clades[:n_transition_clades]) # get just top n_transition_clades clades with highest entropy

    # assign state per clade
    state_dict = {}
    for clade in unique_clades:
        if clade != -1.: #  -1 is not a clade, is just all cells that couldn't be assigned to a clade
            if clade in transition_clades:
                state_dict[clade] = "transitional"
            elif frac_c_dict[clade] >= cutoff:
                state_dict[clade] = "cancer"
            elif frac_nc_dict[clade] >= cutoff:
                state_dict[clade] = "normal"
            else:
                state_dict[clade] = "unsure"
        else:
            state_dict[clade] = "unassigned" # these cells have not been assigned to a clade

    return state_dict


def get_preliminary_names(tree):
    # name the root
    tree.root.name = "root"
    # give arbitrary unique names to non leaf nodes (leaf nods should have unique names already (eg cell ID))
    i = 0
    for node in tree.get_nonterminals():
        if node.name != "root":
            node.name = f"node_{i}"
        i += 1
    # assert there are no duplicate names in nonterminals
    assert len(set(tree.get_nonterminals())) == len(tree.get_nonterminals())

def main():

    # import adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)

    # build neighbor joining tree from cnv matrix
    # print statements inside function (nj tree, distance matrix)
    tree = build_tree(adata, cnv_score_matrix, distance_metric)

    # Save tree
    print("saving tree to temp file")
    tree.write(os.path.join(output_dir, f"cnv_tree_{os.path.basename(input_data_file).removesuffix('.h5ad')}.nwk"), format="newick")
    
    # change format to biopython (easier to manipulate here)
    print("converting scipy tree to biopython tree")
    tree = Phylo.read(os.path.join(output_dir, f"cnv_tree_{os.path.basename(input_data_file).removesuffix('.h5ad')}.nwk"), format="newick")

    # name tree nodes and assert uniqueness (downstream functions depend on unqiue names for all nodes)
    print("assigning unqiue names to tree nodes")
    get_preliminary_names(tree)

    # annotate adata with clades
    print("annotating adata with clades")
    clade_list = divide_tree(tree, n_clades) # list of string tuples, each tuple contains cell ids that belong to same clade
    fixed_clade_list = fix_clades(tree, clade_list) # there might be some clades the surround other clades (i.e. the have "gaps"), this function splits those clades at the gaps
    # adata.obs["cnv_clade"] = None  !!! DO NOT INITIALIZE WITH NONE, this will lead to issues with saving to h5ad!!!
    for i, clade in enumerate(fixed_clade_list):
        print(f"annotating clade {i}")
        mask = adata.obs_names.isin(clade)
        adata.obs.loc[mask, "cnv_clade"] = int(i) # values 0 to n_clades - 1
    adata.obs["cnv_clade"] = adata.obs["cnv_clade"].fillna(-1) # replace cells with no clade by clade = -1 (so downstream throws no errors)

    # annotate adata with cancer_state_inferred_tree (normal, transitional, cancer)
    print("annotating adata with cancer_state_inferred_tree")
    states_dict = get_states(adata, metric, n_transition_clades, cutoff)
    for clade, state in states_dict.items():
        mask = adata.obs["cnv_clade"] == clade
        adata.obs.loc[mask, "cancer_state_inferred_tree"] = state

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + output_dir)


if __name__ == "__main__":
    
    input_data_file, output_dir, cnv_score_matrix, distance_metric, n_clades, metric, n_transition_clades, cutoff, verbose = hf.import_cmd_args(8)
    vprint = hf.make_vprint(verbose)

    main()