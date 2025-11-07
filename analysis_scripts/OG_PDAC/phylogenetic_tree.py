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

def build_tree(adata, cnv_score_matrix):# get cnv score matrix
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
    D = squareform(pdist(X_jittered, metric='correlation')) # correlation = how similar cnv profiles are, invariant to magnitude of expression
    vprint(f"Percentage of NaNs in distance matrix: {np.isnan(D).sum() / D.size * 100}%")
    vprint(f"Number of NaNs in distance matrix: {np.isnan(D).sum()}")
    dm = DistanceMatrix(D, ids=labels)

    # Build neighbor-joining tree
    print("building neighbor-joining tree")
    tree = nj(dm)
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
                clade_leaf_sets.append(tuple(leaves)) # tuple of current clades leaves
                for leaf in clade_leaf_sets[-1]: # newest entry in clade list = this entry
                    leaf_pool.remove(leaf)
                clade_root = node # set nonlocal variable "clade_root" to chosen node
                return True # true means this node is chosen
        
        parent_chosen = recurse(parents[node.name], n_leaves) # recursive call on parent

        if parent_chosen == False: # if the parent was not chosen (because this node is closer to ideal clade size)
            clade_leaf_sets.append(tuple(leaves))
            for leaf in clade_leaf_sets[-1]:
                leaf_pool.remove(leaf)
            clade_root = node # set nonlocal clade root to chosen node
        return None # None means this node OR one of its ancestors was chosen
    

    # sanity check
    if n_clades < 1:
        raise ValueError("n_clades must be at least 1")
    if n_clades == 1:
        return [tuple(tip.name for tip in tree.tips())]
    
    # local vars
    terminals = list(tree.get_terminals())
    parents = all_parents(tree)
    leaf_pool = list(tree.get_terminals(order = "postorder")) # list of available leaves (names of the cells) for selection
    n_leaves = len(leaf_pool) # total number of leaves in tree
    target_clade_size = n_leaves / n_clades # target number of leaves in ideal clade
    clade_leaf_sets = [] # list of chosen clades (list of str tuples)

    # find start node
    i = 0
    while len(clade_leaf_sets) < n_clades and len(leaf_pool) > 1: # repeat until n_clades reached or no leaves left
        print(f"Building clade {str(i)}")
        visitation_order = get_visitation_order_width(tree)
        start_node = find_start_node(tree, visitation_order) # find start node
        clade_root = None # root of clade (will be set in recurse)
        recurse(start_node, 0) # recurse upward from that node and add entry to clade_leaf_sets
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


def get_states(adata, metric, n_transition_clades, cutoff): 
    # metric is one of adata.obs["cancer_state"] or adata.obs["cancer_state_inferred"]
    # cutoff is eg: calde has 90% c 10% nc -> classified as c, 89% c 11% nc, classified as unsure 
    
    # define entropy dict
    # for each clade in adata cnv_clades
        # get the number of c and nc cells according to metric
        # compute entropy
        # add entropy to dict (key = clade, value = entropy)
    # define state dict
    # choose n_transition_clades highest entropy clades and set the value of those clades in state dict to "transitional"
    # set value of remaining clades with > cutoff nc cells as normal
    # set value of remaining clades with > cutoff c cells as cancer
    # return state dict

    clades = np.asarray(adata.obs["cnv_clade"])
    metrics = np.asarray(adata.obs[metric])

    unique_clades = np.unique(clades)

    entropy_dict = {}
    frac_c_dict = {}
    frac_nc_dict = {}

    for clade in unique_clades:
        mask = clades == clade
        clade_metrics = metrics[mask]

        n_total = len(clade_metrics)
        n_c = np.sum(clade_metrics == "c")
        n_nc = np.sum(clade_metrics == "nc")

        if n_total == 0:
            frac_c = 0.0
            frac_nc = 0.0
            entropy = 0.0
        else:
            frac_c = n_c / n_total
            frac_nc = n_nc / n_total
            # Shannon entropy
            if frac_c in (0.0, 1.0):
                entropy = 0.0
            else:
                entropy = -(frac_c * np.log2(frac_c) + frac_nc * np.log2(frac_nc))

        entropy_dict[clade] = entropy
        frac_c_dict[clade] = frac_c
        frac_nc_dict[clade] = frac_nc

    # choose clades with highest entropy
    sorted_clades = sorted(entropy_dict, key=entropy_dict.get, reverse=True)
    transition_clades = set(sorted_clades[:n_transition_clades])

    # assign state per clade
    state_dict = {}
    for clade in unique_clades:
        if clade in transition_clades:
            state_dict[clade] = "transitional"
        elif frac_c_dict[clade] >= cutoff:
            state_dict[clade] = "cancer"
        elif frac_nc_dict[clade] >= cutoff:
            state_dict[clade] = "normal"
        else:
            state_dict[clade] = "unsure"

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

    # Debugging
    """# build neighbor joining tree from cnv matrix
    tree = build_tree(adata, cnv_score_matrix)

    # Save tree
    print("saving tree to temp file")
    tree.write(os.path.join(output_dir, f"cnv_tree{os.path.basename(input_data_file).removesuffix('.h5ad')}.nwk"), format="newick")"""
    
    # change format to biopython (easier to manipulate here)
    tree = Phylo.read(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\PDAC_ductal_cnv_tree.nwk", format="newick")

    # name tree nodes
    get_preliminary_names(tree)

    # annotate adata with clades
    print("annotating adata with clades")
    clade_list = divide_tree(tree, n_clades) # list of string tuples, each tuple contains cell ids that belong to same clade
    adata.obs["cnv_clade"] = None
    for i, clade in enumerate(clade_list):
        mask = adata.obs_names.isin(clade)
        adata.obs.loc[mask, "cnv_clade"] = i # values 0 to n_clades - 1


    # annotate adata with cancer_state_inferred_tree (normal, transitional, cancer)
    print("annotating adata with cancer_state_inferred_tree")
    states_dict = get_states(adata, metric, n_transition_clades, cutoff)
    adata.obs["cancer_state_inferred_tree"] = None
    for clade, state in states_dict.items():
        mask = adata.obs["clade"] == clade
        adata.obs.loc[mask, "cancer_state_inferred_tree"] = state

    # save adata
    print("sacing adata to temp file")
    adata.write(os.path.join(output_dir, os.path.basename(input_data_file)))
    
    # notify subrocess executor of toutput (copies everything in this folder in case of permanent storage)
    print("Output: " + output_dir)


if __name__ == "__main__":
    
    input_data_file, output_dir, cnv_score_matrix, n_clades, metric, n_transition_clades, cutoff, verbose = hf.import_cmd_args(8)
    vprint = hf.make_vprint(verbose)

    main()