# generate a phylogenetic tree from X_CNV (so the cnvs are the evolutionary steps)
# show the tree
# save in format where we can select clades corresponding to cells in adata so we can define transition cells

import scanpy as sc
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj
import helper_functions as hf
import seaborn as sns
import os
import plotly.graph_objects as go
from Bio import Phylo
import json
import re

def visualize_tree(adata):

    # ==========================
    # Configuration
    # ==========================

    tree_file = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\test trees/1500TNtree.nwk"      # output of scikit-bio neighbor-joining
    metadata_columns = obs_columns  # columns from adata.obs to show as rings
    ring_spacing = 0.025            # distance between concentric rings
    base_radius = 1.0               # radius of leaf circle
    print(metadata_columns)

    # ==========================
    # Step 1: Load tree + metadata
    # ==========================

    print("[INFO] Loading tree and metadata...")
    tree = Phylo.read(tree_file, "newick")
    obs = adata.obs.copy()

    print(type(tree))

    # debug
    # get_max_tree depth
    def max_depth(node):
        if not node.clades:  # leaf
            return 0
        return 1 + max(max_depth(c) for c in node.clades)

    maximum_depth = max_depth(tree.root)
    print("Max depth (edges from root):", maximum_depth)


    # ==========================
    # Step 2: Get ordered leaves (no overlaps)
    # ==========================

    def ordered_terminals(tree):
        """Return leaf nodes in topologically consistent postorder."""
        return list(tree.get_terminals(order="postorder"))

    terminals = ordered_terminals(tree)
    n_leaves = len(terminals)
    leaf_names = [leaf.name for leaf in terminals]

    print(f"[INFO] Found {n_leaves} leaves (cells).")
    print(f"dytpe of terminals: {type(terminals)}")
    print(f"dtype of elemts of terminals: {type(terminals[0])}")
    print(f"[INFO] found terminal leaves: {terminals}")


    # ==========================
    # Step 3: Assign radial coordinates
    # ==========================

    # give arbitrary unique names to non leaf nodes
    i = 0
    for node in tree.get_nonterminals():
        node.name = f"node_{i}"
        i += 1
    # assert there are no duplicate names in nonterminals
    assert len(set(tree.get_nonterminals())) == len(tree.get_nonterminals())

    def all_parents(tree):
        parents = {}
        for clade in tree.find_clades(order="level"):
            for child in clade:
                parents[child.name] = clade
        return parents
    parents = all_parents(tree) # key = node name, value = its parent's node name

    # output as json
    """with open("parents.json", "w") as f:
        json.dump(parents, f, indent=4)"""
        

    def compute_radial_layout(tree, leaf_names, terminal_node_radius):
        """Compute (theta, r) coordinates for all nodes in the tree."""
        coords = {}
        for terminal in terminals:
            coords[terminal.name] = (0, terminal_node_radius)

        # recursively assign internal node radii inward based on depth
        def set_internal(node, depth, angle): # angle is the total sum of all angles up the current node on the entire tree (depth first)
            if node.name not in leaf_names: # so all terminal nodes will be skipped (so we need not check again below)
                if node.name != "root":
                    
                    # NAME THE NODE
                    # define local vars
                    parent = parents[node.name] # get the node's parent (clade object)
                    siblings = list(parent.clades) # get the node's parent's children (i.e. siblings + self) (clade objects)
                    siblings.remove(node) # now remove the node itself
                    pattern = rf"^{re.escape(str(parent.name))}_\d+$" # this is what a proper sibling name looks like f"parent_{i}"
                    # count children of the node's parent with proper names
                    properly_named_siblings = 0
                    for element in siblings:
                        if re.match(pattern, element.name):
                            properly_named_siblings += 1 
                    # assign a proper name to the node
                    node.name = f"{str(parents[node.name].name)}_{properly_named_siblings + 1}" # name the node as the nth child of its parent (n = properly_named_siblings + 1, so if a parent has 5 childern, then it will be parent_1 for the first one that is checked, then parent_2, and so on up to root_1_1_1_2_5_1_7_4_2_2...)
                
                # ANGLE THE NODE
                # angle the current node itself
                coords[node.name] = (angle, (base_radius/maximum_depth) * depth)
                
                # CHECK FOR TERMINAL CHILDREN
                # check for existence
                terminal_children = list(set(node.clades).intersection(set(terminals))) # children of node that are leaf nodes
                # angle the terminal children if they exist
                added_angle = 0 # define here to avoid unbound variable error
                if terminal_children != []:
                    angled_terminal_children = 0
                    for tc in terminal_children:
                        _, og_radius = coords[tc.name] # get tuple values
                        coords[tc.name] = (angle + ((2*np.pi)/n_leaves)*angled_terminal_children, og_radius) # angle = to what the current node itself gets + modifer per terminal child of current node
                        angled_terminal_children += 1
                    # determine angle added because of terminal children
                    added_angle = angled_terminal_children * ((2*np.pi)/n_leaves) # this is the total angle used by the terminal children, it should be added to that of each proper child to avoid overlaps

                # TOTAL ANGLE
                total_angle = angle + added_angle # combination of prior angle + gained angle from terminal children

                # RECURSIVE CALL
                for child in node.clades:
                    if child not in terminal_children:
                        total_angle = set_internal(child, depth + 1, total_angle)

                # RETURN (this can only trigger once all recurscive calls are done)
                return total_angle
                    
                
        tree.root.name = "root"
        coords["root"] = (0, 0)
        set_internal(tree.root, 0, 0)

        return coords

    coords = compute_radial_layout(tree, leaf_names, base_radius) # changes names, so parents dict must be recomputed


    def lower_branch_node_radii(tree, node):

        # define local consts
        children = list(node.clades)
        non_terminal_children = list(set(children).difference(set(terminals)))

        # recurse on all non terminal children
        child_radii = [] # list to store radii of all non terminal children
        if non_terminal_children != []:
            for ntc in non_terminal_children:
                child_radius = lower_branch_node_radii(tree, ntc)
                child_radii.append(child_radius)
        min_child_radius = min(child_radii) if child_radii != [] else base_radius # if there are only terminal children, the will all have the base radius

        if node.name != "root": # avoid moving the root node
            # move node (all nodes descending from it have already been moved)
            og_angle, og_radius = coords[node.name]
            new_radius = min_child_radius - base_radius/maximum_depth # so say the child that is closest to the center is at 1 and 50 nodes then this nodes at 1 - 1/50 = 49/50
            coords[node.name] = (og_angle, new_radius)

            return new_radius        

        
    # move internal nodes closer to terminal nodes for visibility
    lower_branch_node_radii(tree, tree.root)


    #DEBUGGING: print tree structure to json
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child.name] = str(clade)
    with open("parents_proper_names.json", "w") as f:
        json.dump(parents, f, indent=4)

    # convert polar → cartesian
    xy_coords = {}
    for name, (theta, r) in coords.items():
        xy_coords[name] = (r * np.cos(theta), r * np.sin(theta))


    # ==========================
    # Step 4: Define color mappings for metadata
    # ==========================

    def map_colors(df, column):
        cats = df[column].unique()
        palette = sns.color_palette("tab10", len(cats)).as_hex()
        cmap = dict(zip(cats, palette))
        return df[column].map(cmap), cmap

    column_color_maps = {}
    for col in metadata_columns:
        obs[f"{col}_color"], cmap = map_colors(obs, col)
        column_color_maps[col] = cmap


    # ==========================
    # Step 5: Draw edges (gray lines)
    # ==========================

    edge_traces = []
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            x0, y0 = xy_coords.get(clade.name, (0, 0))
            x1, y1 = xy_coords.get(child.name, (0, 0))
            edge_traces.append(
                go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color="lightgray", width=1),
                    hoverinfo="none",
                    showlegend=False,
                )
            )


    # ==========================
    # Step 6: Draw colored annotation rings
    # ==========================

    base_coords = np.array([xy_coords[n] for n in leaf_names])
    node_traces = []

    for i, col in enumerate(metadata_columns):
        colors = ["black" for _ in leaf_names]                                  # obs.loc[leaf_names, f"{col}_color"]
        xs, ys = base_coords[:, 0], base_coords[:, 1]
        scale = 1 + ring_spacing * i  # expand outward for each ring

        node_traces.append(
            go.Scatter(
                x=xs * scale,
                y=ys * scale,
                mode="markers",
                marker=dict(size=4, color=colors, line=dict(width=0)),
                name=col,
                hovertext=[
                  f"{leaf}"                        #f"{leaf}: {col}={obs.loc[leaf, col]}"
                  for leaf in leaf_names                           # for leaf in leaf_names
                ],
                hoverinfo="text",
            )
        )

    # DEBUG
    # write coords and xycoords to json
    with open("coords.json", "w") as f:
        json.dump(coords, f, indent=3)
    with open("xycoords.json", "w") as f:
        json.dump(xy_coords, f, indent=3)
    print("coords written to json")

    # draw point at (0,0) for debugging
    node_traces = []
    node_traces.append(
        go.Scatter(
            x=[0],
            y=[0],
            mode="markers",
            marker=dict(size=20, color="red", line=dict(width=0)),
        )
    )

    # draw point at tree root location
    print(f"Tree root is at:\n{coords[tree.root.name]}")
    node_traces.append(
        go.Scatter(
            x= [int(xy_coords[tree.root.name][0])],   
            y= [int(xy_coords[tree.root.name][1])],
            mode="markers",
            marker=dict(size=20, color="blue", line=dict(width=0)),
        )
    )

    # draw smaller point at each internal node
    for clade in tree.find_clades(order="level"):
        x, y = xy_coords.get(clade.name, (0, 0))
        node_traces.append(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers",
                marker=dict(size=5, color="black", line=dict(width=0)),
                hovertext=[clade.name],
                hoverinfo="text",
            )
        )

    # ==========================
    # Step 7: Combine and render
    # ==========================

    fig = go.Figure(edge_traces + node_traces)  # ADD + node_traces after debugging

    fig.update_layout(
        showlegend=True,
        title="CNV Tree with Multi-Ring Annotations",
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        plot_bgcolor="white",
        width=900,
        height=900,
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    fig.show()

    # optional: save to file
    # fig.write_image("cnv_tree_multiring.svg")
    

def main():

    # import adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)

    # get cnv score matrix
    vprint("getting cnv score matrix")
    X = adata.obsm[cnv_score_matrix].toarray() # need to turn scipy sparse csr matrix into numpy array for pdist to work

    # DEBUGGING, print percentage of NaNs in X
    print(f"Percentage of NaNs in X: {np.isnan(X).sum() / X.size * 100}%")
    print(f"Number of NaNs in X: {np.isnan(X).sum()}")
    print(f"X: {X}")

    labels = adata.obs_names.tolist()

    # add a tiny jitter so pearson correlation can be computed for all cells 
    # (if there are rows with 0 variance (ie all 0s) then pearson correlation cannot be computed --> NaNs)
    epsilon = 1e-16
    vprint("adding jitter to X")
    X_jittered = X + epsilon * np.random.randn(*X.shape)

    # Compute distance matrix
    print("computing distance matrix")
    D = squareform(pdist(X_jittered, metric='correlation')) # correlation = how similar cnv profiles are, invariant to magnitude of expression
    print(f"Percentage of NaNs in X: {np.isnan(D).sum() / D.size * 100}%")
    print(f"Number of NaNs in X: {np.isnan(D).sum()}")
    print(f"X: {D}")
    dm = DistanceMatrix(D, ids=labels)

    # Build neighbor-joining tree
    print("building neighbor-joining tree")
    tree = nj(dm)

    # Save tree
    print("saving tree to temp file")
    tree.write(os.path.join(output_dir, f"cnv_tree{os.path.basename(input_data_file)}.nwk"), format="newick")
    tree.write(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\test.nwk", format="newick")

    # Visualize tree
    print("visualizing tree")
    visualize_tree(adata)

if __name__ == "__main__":
    
    # input_data_file, output_dir, cnv_score_matrix, obs_columns, show, verbose = hf.import_cmd_args(5)
    # vprint = hf.make_vprint(verbose)

    input_data_file = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\reduced\reduced_PDAC_ductal_cell.h5ad"
    output_dir = "None"
    cnv_score_matrix = None
    obs_columns = ["cancer_state"]
    show = True
    verbose = True
    vprint = hf.make_vprint(verbose)

    # debugging
    print(obs_columns)

    import sys
    adata = sc.read_h5ad(input_data_file)
    visualize_tree(adata)
    sys.exit(0)

    main()