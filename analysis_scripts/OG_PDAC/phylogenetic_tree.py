# generate a phylogenetic tree from X_CNV (so the cnvs are the evolutionary steps)
# show the tree
# save in format where we can select clades corresponding to cells in adata so we can define transition cells

import scanpy as sc
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj
from skbio import TreeNode
import helper_functions as hf
import os

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

def divide_tree(tree, n_clades) -> list[tuple[str, ...]]: 

    # iterate through tree from root to children
        # 2 lists, one with current depth nodes, one with their parents
        # check list length (save in dict (key = depth, value = list length))
        # if current depth list length crosses n_clades
            # take the clades in the list (current / parents) where abs(len - n_clades) is lower
            # get leaf nodes in those clades, one tuple per clade, add tuples to list
            # return list

    # sanity check
    if n_clades < 1:
        raise ValueError("n_clades must be at least 1")
    if n_clades == 1:
        return [tuple(tip.name for tip in tree.tips())]

    # level-order traversal
    current_level = [tree]
    parent_level = []
    depth = 0
    level_counts = {}

    # Traverse levels until we exceed or reach n_clades
    while current_level:
        level_counts[depth] = len(current_level)
        if len(current_level) >= n_clades:
            break
        # prepare next level
        parent_level = current_level
        next_level = []
        for node in current_level:
            next_level.extend(node.children)
        current_level = next_level
        depth += 1

    # choose between this level and the previous one, whichever is closer
    if not parent_level:
        best_level = current_level
    else:
        if abs(len(current_level) - n_clades) < abs(len(parent_level) - n_clades):
            best_level = current_level
        else:
            best_level = parent_level

    # gather all leaves under each node at the best level
    clade_leaf_sets = []
    for node in best_level:
        leaves = tuple(tip.name for tip in node.tips())
        clade_leaf_sets.append(leaves)

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
    
    tree = TreeNode.read(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\PDAC_ductal_cnv_tree.nwk")

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
    
    input_data_file, output_dir, cnv_score_matrix, n_clades, metric, n_transition_clades, cutoff, verbose = hf.import_cmd_args(6)
    vprint = hf.make_vprint(verbose)

    main()