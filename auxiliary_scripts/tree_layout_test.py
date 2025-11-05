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
from typing import Literal

def visualize_tree(adata, tree_file, obs_columns: list = [], ring_spacing = 0.025, base_radius = 1, target_circumference = 2*np.pi * 0.9, sort_order: Literal["lowest_width_first", "lowest_depth_first"] = "lowest_width_first", show = True, save = False, debug = False):

    # PRELIMINARY NAMING OF NODES
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


    # GET ORDERED LEAVES
    def ordered_terminals(tree):
        """Return leaf nodes in topologically consistent postorder."""
        return list(tree.get_terminals(order="postorder"))
    

    # FIND MAX DEPTH (needed for computing node radii in radial layout)
    def max_depth(node):
        if not node.clades:  # leaf
            return 0
        return 1 + max(max_depth(c) for c in node.clades)
    

    # GET VISITATION ORDER BY CLADE DEPTH
    # This is a dictionary of tuples
    # the key is a node's name, and the value is a list of tuples
    # each tuple corresponds to one of the node's children
    # the first entry in the tuple is the child (clade object)
    # the second entry is the distance of that child to its closest descendant that only has terminal children
    # i.e. if a node has a child that only has terminal children, then the second entry in the tuple 
    # if a node only has terminal children, then the second entry in the tuple is 0
    # if a node has a child that only has terminal children, then the second entry in the tuple is 1
    # and so on
    # The list is sorted in ascending order by distance
    # this means that children with the lowest such distance are first in the list (and will be visited first in the radial layout
    # (this leads to the radial layout mostly avoiding overlapping branches)
    # IMPORTANT: this is invalidated if the nodes are renamed
    def get_visitation_order_depth(tree):
        visitation_order = {}

        def recurse(node):
            children = list(node.clades)
            non_terminal_children = list(set(children).difference(set(terminals)))
            if non_terminal_children == []:
                return 0 # if the node itself only has terminaly children, its distance to such a nod is evidently 0 (bcs it is such a node)

            distance_list = []
            for ntc in non_terminal_children:
                distance_to_last_ntc = recurse(ntc) # this is the distance from the ntc to its closest descendant that only has terminal children (0 = the node itself only has terminal children, 1 = a directly connected node only has terminal children, etc.)
                distance_list.append((ntc, distance_to_last_ntc)) # tuple of child as clade object and its above described distance

            distance_list.sort(key=lambda x: x[1]) # sort in ascending order by distance
            visitation_order[node.name] = distance_list

            return distance_list[0][1] + 1
        
        recurse(tree.root) # fills visitation_order dict

        return visitation_order
    
    # ALTERNATIVE: GET VISITATION ORDER BY CLADE WIDTH
    # wider clade as are visited liast
    # more terminal nodes in clade mines it is wider
    # also returns a dict of list(tuple, tuple, ...)
    # almost everything is identical to get_visitation_order_depth
    # but the number in the tuple refers to the total number of terminal nods in the child clade
    # i.e. if a node (key) has a child clade that only has 1 terminal child, then the second entry in the tuple is 1
    # tuples are sorted in ascending order by width for each key
    def get_visitation_order_width(tree):
        visitation_order = {}
        for clade in tree.find_clades(order="level"):
            children = list(clade.clades)
            non_terminal_children = list(set(children).difference(set(terminals)))

            width_list = []
            for ntc in non_terminal_children:
                width = len(list(ntc.get_terminals())) # width of the clade
                width_list.append((ntc, width)) # tuple of child as clade object and its width
            
            width_list.sort(key=lambda x: x[1]) # sort in ascending order by width
            visitation_order[clade.name] = width_list
        
        return visitation_order


    # GET PARENTS DICTIONARY
    # the keys are node names, and the values are the clade object corresponding to its parent
    # IMPORTANT: this is invalidated if nodes are renamed (or maybe not, but radial layout function works as is)
    def all_parents(tree):
        parents = {}
        for clade in tree.find_clades(order="level"):
            for child in clade:
                parents[child.name] = clade
        return parents
    

    # COMPUTE RADIAL LAYOUT OF TREE
    # assigns polar coordinates to ever node in the tree to create a circular tree structure for plotting
    # recursive function -> avoid heavy computation within (this is why we precompute parents and visitation orders -> retrieve before renaming of node to stay valid)
    def compute_radial_layout(tree, leaf_names, terminal_node_radius):
        """Compute (theta, r) coordinates for all nodes in the tree."""
        coords = {}
        for terminal in terminals:
            coords[terminal.name] = (0, terminal_node_radius) # set radius for all terminal nodes, angle is set later

        # recursively assign internal node radii and angles based on related nodes (should start from root)
        def set_internal(node, depth, angle): # angle is the total sum of all angles up the current node on the entire tree (so we avoid all collisions with previous nodes)
            if node.name not in leaf_names: # so all terminal nodes will be skipped (so we need not check again below)
                
                # before renaming, extract visitation order (remaning breaks the mapping)
                if any(child not in terminals for child in node.clades): # if the current node only has terminal children, it will not have an entry in visitation_order
                    # extract visitation order for direct children of the current node
                    local_visitation_order = visitation_order[node.name] # list of tuples, lower index = visit first (child (clade object), distance as described in get_visitation_order)

                if node.name != "root": # we do not want to rename the root node
                    
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
                        coords[tc.name] = (angle + ((target_circumference)/n_leaves)*angled_terminal_children, og_radius) # angle = to what the current node itself gets + modifer per terminal child of current node (assign new tuple as tuples are immutable)
                        angled_terminal_children += 1
                    # determine angle added because of terminal children
                    added_angle = angled_terminal_children * ((target_circumference)/n_leaves) # this is the total angle used by the terminal children, it should be added to that of each proper child to avoid overlaps

                # TOTAL ANGLE
                total_angle = angle + added_angle # combination of prior angle + gained angle from terminal children

                # RECURSIVE CALL
                if any(child not in terminals for child in node.clades): # need to make sure that the node has non terminal children (visitation order contains no keys for terminal nodes)
                    for child, distance in local_visitation_order: # already ordered so children with smaller distance to closest descendant with only terminal children are visited first
                        if child not in terminal_children: # do not visit terminal nodes / leaf nodes
                            total_angle = set_internal(child, depth + 1, total_angle) # total angle increased every time a a node with terminal children is visited

                # RETURN (this can only trigger once all recurscive calls are done for the current node)
                return total_angle
                    
        coords["root"] = (0, 0) # seet base coords for root
        set_internal(tree.root, 0, 0) # call the function on the root
        return coords
    

    # MOVE INTERNAL NODES CLOSER TO TERMINAL NODES FOR VISIBILITY
    def increase_branch_node_radii(tree, node):

        # define local consts
        children = list(node.clades)
        non_terminal_children = list(set(children).difference(set(terminals)))

        # recurse on all non terminal children
        child_radii = [] # list to store radii of all non terminal children
        if non_terminal_children != []:
            for ntc in non_terminal_children:
                child_radius = increase_branch_node_radii(tree, ntc)
                child_radii.append(child_radius)
        min_child_radius = min(child_radii) if child_radii != [] else base_radius # if there are only terminal children, they will all have the base radius

        if node.name != "root": # avoid moving the root node
            # move node (all nodes descending from it have already been moved)
            og_angle, og_radius = coords[node.name]
            new_radius = min_child_radius - base_radius/maximum_depth # so say the child that is closest to the center is at 1 and 50 nodes then this nodes at 1 - 1/50 = 49/50
            coords[node.name] = (og_angle, new_radius)

            return new_radius   


    # DEFINE COLOR MAPS FOR OBS COLUMNS
    def map_colors(df, column):
        cats = df[column].unique()
        palette = sns.color_palette("tab10", len(cats)).as_hex()
        cmap = dict(zip(cats, palette))
        return df[column].map(cmap), cmap
    

    # DRAW EDGES AS GREY LINES
    # draws an edge between any 2 connected nodes
    def add_edges(tree, edge_traces, xy_coords):
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
    

    # DRAW ANNOTATION RINGS
    # draws one ring per obs column
    # different categories get different colors
    def add_annotation_rings(node_traces, obs_columns, xy_coords, leaf_names):
        base_coords = np.array([xy_coords[n] for n in leaf_names])
        if obs_columns != []:
            vprint(f"Obs columns to draw annotation rings for: {obs_columns}")
            for i, col in enumerate(obs_columns):
                colors = obs.loc[leaf_names, f"{col}_color"]
                xs, ys = base_coords[:, 0], base_coords[:, 1]
                scale = 1 + ring_spacing * i  # expand outward for each ring

                vprint(f"Drawing annotation ring for {col}")
                node_traces.append(
                    go.Scatter(
                        x=xs * scale,
                        y=ys * scale,
                        mode="markers",
                        marker=dict(size=4, color=colors, line=dict(width=0)),
                        name=col,
                        hovertext=[
                        f"{leaf}, {col}={obs.loc[leaf, col]}"
                        for leaf in leaf_names                           
                        ],
                        hoverinfo="text",
                    )
                )

    # LOAD TREE AND METADATA
    vprint("Loading tree and obs columns...")
    tree = Phylo.read(tree_file, "newick")
    obs = adata.obs.copy()

    # ASSIGN PRELIMINARY NAMES (needed for functions that work on the tree)
    vprint("Assigning preliminary names to tree nodes...")
    get_preliminary_names(tree)


    # GET MORE GENERAL TREE PARAMTERS (needed for function behavior)
    vprint("Finding leaves...")
    terminals = ordered_terminals(tree)
    n_leaves = len(terminals)
    leaf_names = [leaf.name for leaf in terminals]
    vprint(f"Found {n_leaves} leaves (cells) in tree")
    # GENERAL TREE PARAMTER
    vprint("Finding maximum depth of tree...")
    maximum_depth = max_depth(tree.root)
    vprint("Max depth of tree (n_nodes from root):", maximum_depth)

    # VISITATION ORDER (needed for tree traversal and formatting)
    vprint("Computing visitation order dictionary for tree traversal...")
    if sort_order == "lowest_depth_first":
        visitation_order = get_visitation_order_depth(tree) # use for set internal in compute_radial_layout
    elif sort_order == "lowest_width_first":
        visitation_order = get_visitation_order_width(tree)

    # PARENT CHILD RELATIONSHIPS (needed for function behavior, do not move this, names are changed in compute_radial_layout)
    vprint("Creating dictionary of child to parent mappings...")
    parents = all_parents(tree) # key = node name, value = its parent's node name

    # GET ACTUAL TREE LAYOUT (node locations, uses recursion to traverse the tree according to visitation order)
    vprint("Computing radial layout of tree nodes...")
    coords = compute_radial_layout(tree, leaf_names, base_radius) # changes names, so parents dict must be recomputed

    # INCREASE BRANCH NODE RADIUS (usually the are quite centered, which makes branches quite dense and hard to see)
    vprint("Increasing branch node radii for visibility...")
    increase_branch_node_radii(tree, tree.root)

    # MAPPING FOR ANNOTATION RINGS
    vprint("Mapping obs columns to colors...")
    column_color_maps = {}
    if obs_columns != []:
        for col in obs_columns:
            obs[f"{col}_color"], cmap = map_colors(obs, col)
            column_color_maps[col] = cmap

    # CONVERT POLAR TO XY COORDS
    vprint("Converting polar coordinates to xy coordinates for plotting...")
    xy_coords = {}
    for name, (theta, r) in coords.items():
        xy_coords[name] = (r * np.cos(theta), r * np.sin(theta))

    # DRAW EDGES
    vprint("Drawing edges between connected nodes...")
    edge_traces = []
    add_edges(tree, edge_traces, xy_coords)

    # DRAW ANNOTATION RINGS (1 ring per obs column)
    vprint("Drawing annotation rings...")
    node_traces = []
    add_annotation_rings(node_traces, obs_columns, xy_coords, leaf_names)

    # DRAW POINT AT ROOT LOCATION (for orientation)
    vprint("Drawing point at root location...")
    vprint(f"Tree root is at:\n{coords[tree.root.name]}")
    node_traces.append(
        go.Scatter(
            x= [int(xy_coords[tree.root.name][0])],   
            y= [int(xy_coords[tree.root.name][1])],
            mode="markers",
            marker=dict(size=10, color="black", line=dict(width=0)),
            hovertext=[tree.root.name],
            hoverinfo="text",
        )
    )


    # FIGURE RELATED DEBUGGING POINTS (must be done before rendering as they are drawn on the figure)
    if debug:
        # draw point at (0,0) for debugging
        print("DEBUG: Drawing point at (0,0)")
        node_traces = []
        node_traces.append(
            go.Scatter(
                x=[0],
                y=[0],
                mode="markers",
                marker=dict(size=20, color="red", line=dict(width=0)),
            )
        )

        # draw smaller point at each internal node
        print("DEBUG: Drawing points at each internal node")
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

   
    # RENDER THE FIGURE
    vprint("Rendering figure...")
    fig = go.Figure(edge_traces + node_traces)  # combine edge and node traces
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

    # SHOW AND SAVE FIGURE
    if show: 
        fig.show()
    if save: 
        print("Oh, you want save to save the figure? That's great! But you will have to do it manually, bcs plotly is a bit weird.")
        print(f"You may want to save to this path: {(os.path.join(output_dir, f"{os.path.basename(tree_file).removesuffix('.nwk')}.png"))}")

        if show != True:
            print("I shall now show you the figure, so you may save it.")
            fig.show()
        if show == True: # don't need to show again if it is already shown
            print("You already have the figure shown, so I shan't overtax your RAM by showing it again.")

        
    # FILE RELATED DEBUGGING OPTIONS
    if debug:
        #DEBUGGING: write all parent-child relationships to json (with latest names)
        print("DEBUG: Writing all parent-child relationships to json...")
        parents = {}
        for clade in tree.find_clades(order="level"):
            for child in clade:
                parents[child.name] = str(clade)
        with open("parents_proper_names.json", "w") as f:
            json.dump(parents, f, indent=4)

        # DEBUGGING: write visitation order to json (with lates names)
        print("DEBUG: Writing visitation order to json...")
        if sort_order == "lowest_depth_first":
            visitation_order = get_visitation_order_depth(tree)
        elif sort_order == "lowest_width_first":
            visitation_order = get_visitation_order_width(tree) # dict of lists of tuples (tuple = (calde, number))
        # change clade to string (clade is not json serializable)
        visitation_order_str = {}
        for item in visitation_order.items():
            key = item[0]
            tuple_list = item[1]
            new_tuple_list = []
            for t in tuple_list:
                new_tuple_list.append((str(t[0]), t[1])) # need to turn clade into string for json
            visitation_order_str[key] = new_tuple_list    
        # dump to json for debugging
        with open("visitation_order.json", "w") as f:
            json.dump(visitation_order_str, f, indent=4)

        # DEBUGGING: write coords to json
        print("DEBUG: Writing radial and xycoords to json...")
        with open("coords.json", "w") as f:
            json.dump(coords, f, indent=3)
        with open("xycoords.json", "w") as f:
            json.dump(xy_coords, f, indent=3)





    

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
    tree.write(os.path.join(output_dir, f"cnv_tree{os.path.basename(input_data_file).removesuffix('.h5ad')}.nwk"), format="newick")
    tree.write(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\test.nwk", format="newick")

    # Visualize tree
    print("visualizing tree")
    visualize_tree(adata)

if __name__ == "__main__":
    
    # input_data_file, output_dir, cnv_score_matrix, obs_columns, show, verbose = hf.import_cmd_args(5)
    # vprint = hf.make_vprint(verbose)

    input_data_file = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\reduced\reduced_PDAC_ductal_cell.h5ad"
    output_dir = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree"
    cnv_score_matrix = None
    obs_columns = ["cancer_state", "cancer_state_inferred"]
    show = True
    verbose = False
    vprint = hf.make_vprint(verbose)

    import sys
    adata = sc.read_h5ad(input_data_file)
    visualize_tree(adata, tree_file=r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\PDAC_ductal_cnv_tree.nwk", save=True, obs_columns=obs_columns)
    sys.exit(0)

    main()