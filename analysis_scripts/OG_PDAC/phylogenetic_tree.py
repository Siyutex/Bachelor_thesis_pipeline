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

    tree_file = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\50TNtree.nwk"      # output of scikit-bio neighbor-joining
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
        

    def compute_radial_layout(tree, leaves):
        """Compute (theta, r) coordinates for all nodes in the tree."""
        coords = {}

        # assign leaves evenly around the circle
        leaf_angles = np.linspace(0, 2 * np.pi, len(leaves), endpoint=False)
        for i, leaf in enumerate(leaves):
            coords[leaf.name] = (leaf_angles[i], base_radius)

        # recursively assign internal node radii inward based on depth
        def set_internal(node, depth, iterator=None):
            if not iterator:
                iterator = 0
            if node.name not in coords: # so all terminal nodes will be skipped (so we need not check again below)
                if node.name != "root":
                    parent = parents[node.name] # get the node's parent (clade object)
                    children = list(parent.clades) # get the node's parent's children (i.e. siblings + self) (clade objects)
                    pattern = rf"^{re.escape(str(parent.name))}_\d+$" # this is what a proper sibling name looks like f"parent_{i}"
                    # count children of the node's parent with proper names
                    properly_named_children = 0
                    for child in children:
                        if re.match(pattern, child.name):
                            properly_named_children += 1
                    node.name = f"{str(parents[node.name].name)}_{properly_named_children + 1}" # name the node as the nth child of its parent (n = properly_named_siblings + 1, so if a parent has 5 childern, then it will be parent_1 for the first one that is checked, then parent_2, and so on up to root_1_1_1_2_5_1_7_4_2_2...)
                

                    # for terminal children of the node's parent, add the angle
                    angled_terminal_children = 0
                    terminal_children = list(set(children).intersection(set(terminals)))
                    for terminal_child in terminal_children:
                        _, og_radius = coords[terminal_child.name] # get tuple values
                        parent_angle = coords[parent.name][0] # this is the angle of the parten
                        coords[terminal_child.name] = (parent_angle + (360/n_leaves)*angled_terminal_children, og_radius)
                        angled_terminal_children += 1

                    added_angle = angled_terminal_children * (360/n_leaves) # this is the total angle used by the terminal children, it should be added to that of each proper child to avoid overlaps


                    # find suitable angle for the node
                    parent_angle = coords[parent.name][0]
                    angle = parent_angle + (360/n_leaves)*properly_named_children + (360/n_leaves)*angled_terminal_children# TODO the constant needs to be expressed as a function of tree structure (if 9000 terminal nodes then we need to space them by 2 * np.pi/9000 (we are using radians here) each to cover the circle, this should also apply to internal nodes to get a nice tree)
                
                    coords[node.name] = (angle, 0.005 * depth) # if the node name == "root" then angle = 0 and depth = 0 (TODO the constant needs to be expressed as a function of tree structure)

                print(f"Node {node.name} is the {iterator}th internal node. It's radius is {coords[node.name][1]} and its theta is {coords[node.name][0]}")
                
                iterator += 1
            for child in node.clades:
                print(f"Child of node {node.name}: {child.name}")
                set_internal(child, depth + 1, iterator)

        tree.root.name = "root"
        coords["root"] = (0, 0)
        set_internal(tree.root, 0)

        # now fix terminal node angles (TODO make it so that these angles are only set once)
        """for node in terminals:
            parent = parents[node.name] # parent as clade object
            parent_angle = coords[parent.name][0] # this is the angle of the parten
            og_angle, og_radius = coords[node.name] # get tuple values
            coords[node.name] = (parent_angle, og_radius)""" # assign new tuple bcs tuples are immutable

        return coords

    coords = compute_radial_layout(tree, terminals)

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
                  "gode"                           #f"{leaf}: {col}={obs.loc[leaf, col]}"
                                             # for leaf in leaf_names
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
        print(f"drawing point at {x}, {y}")
        node_traces.append(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers",
                marker=dict(size=5, color="black", line=dict(width=0)),
            )
        )

    # ==========================
    # Step 7: Combine and render
    # ==========================

    fig = go.Figure(edge_traces)  # ADD + node_traces after debugging

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
    
    input_data_file, output_dir, cnv_score_matrix, obs_columns, show, verbose = hf.import_cmd_args(5)
    vprint = hf.make_vprint(verbose)

    # debugging
    print(obs_columns)

    import sys
    adata = sc.read_h5ad(input_data_file)
    visualize_tree(adata)
    sys.exit(0)

    main()