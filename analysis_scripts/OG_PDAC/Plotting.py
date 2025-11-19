# this script will produce 2D projections af the input file, either PCA or UMAP depending on the projection parameter

# main libs
import scanpy as sc
import numpy as np
import pandas as pd
# native libs
import json
import sys
import os
import re
from typing import Literal
import warnings
# plotting
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
# phylogenetic tree
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj
from Bio import Phylo
# other
import helper_functions as hf


"""def plot_with_external_legend(adata, projection: str, show=True, root_cell_idx: int = None, **kwargs):

    # Prevent Scanpy from immediately showing the plot  
    if projection == "UMAP":
        sc.tl.umap(adata)
        fig = sc.pl.umap(adata,  show=False, **kwargs)
    elif projection == "PCA":
        fig = sc.pl.pca(adata, show=False, **kwargs)
    ax = fig.axes

    if root_cell_idx != None:
        umap_coords = adata.obsm['X_umap']
        plt.scatter(umap_coords[root_cell_idx, 0], umap_coords[root_cell_idx, 1],
            color='red', s=100, edgecolor='black', label='highlighted cell')
        
        plt.show()
        return


    n_cols = 1
    # more than 15 legend entries -> 2 columns
    if len(ax.get_legend_handles_labels()[0]) > 20:
        n_cols = 2
    
    # Move legend outside and shrink font
    ax.legend(
        bbox_to_anchor=(1.05, 1),  # move right of plot
        loc='upper left',
        fontsize='small',
        title=kwargs.get('color', ''),  # use color key as legend title
        ncol=n_cols
    )
    
    # Adjust layout to avoid clipping
    plt.tight_layout()

    # always save to tempdir (which is output_data_dir)
    plt.savefig(os.path.join(output_data_dir, f"{projection}_colored_by_{kwargs.get('color', '')}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))

    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()"""


def limit_cells(adata, isolation_dict):
    # isolation dict should have key = obs columns in adata, values = list of entries in that column to keep
    # in the end this function will produce the intersection of all entries (eg cell_type: ["ductal"], cancer_state:["normal","transitional"] will isolate cells that or ductal and either normal or transitional)
    
    vprint("copying adata")
    internal_adata = adata.copy()

    for obs_column, entry_list in isolation_dict.items():
        vprint(f"limiting {obs_column} to {entry_list}")
        if obs_column in internal_adata.obs.columns and all(entry in internal_adata.obs[obs_column].values for entry in entry_list): # check if column exists and needed values exist in it
            vprint("limiting possible")
            internal_adata = internal_adata[internal_adata.obs[obs_column].isin(entry_list)].copy()
        elif obs_column not in internal_adata.obs.columns:
            vprint("obs column not found")
            warnings.warn(f"obs column {obs_column} not found in adata.obs, skipping to next column...")
        elif not all(entry in internal_adata.obs[obs_column].values for entry in entry_list):
            for entry in entry_list:
                if entry not in internal_adata.obs[obs_column].values:
                    vprint(f"entry {entry} not found in {obs_column}")
                    warnings.warn(f"entry {entry} not found in adata.obs['{obs_column}'], skipping to next column...")

    return internal_adata


def matrix_to_anndata(adata, matrix_key):
    """createas a new anndata with the given matrix (any layer, obsm, or X) as adata.X and returns it.
       preserves obs annotations, only presever var annotatiosn if shape[1] matches
       
       DOES NOT MODIFY ADATA"""

    # check if matrix_key is unique among adata.layers, adata.obsm, and adata.X
    keys = list(adata.layers.keys()) + list(adata.obsm.keys()) + ["X"]
    if keys.count(matrix_key) > 1:
        raise ValueError(f"Matrix key '{matrix_key}' is not unique among layers, obsm, and X")

    # assign correct matrix
    if matrix_key    == "X":
        mat = adata.X.copy()
    elif matrix_key in adata.layers:
        mat = adata.layers[matrix_key].copy()
    elif matrix_key in adata.obsm:
        mat = adata.obsm[matrix_key].copy()
    else:
        raise ValueError(f"Key '{matrix_key}' not found in X, layers, or obsm")
    
    # assign obs and var
    obs = adata.obs.copy()
    var = adata.var.copy()

    # create new anndata, assign var if possible
    if len(var) == mat.shape[1]:
        vprint("Obsm matrix feature count matches adata.var, assigning var and obs to obsm_adata...")
        matrix_adata = sc.AnnData(X=mat, obs=obs, var=var)
    else:
        vprint("Obsm matrix feature count does not match adata.var, only assigning obs to obsm_adata...")
        matrix_adata = sc.AnnData(X=mat, obs=obs)

    return matrix_adata

def get_DEGs(adata, groupby):

    # rank genes groups overrides any previous contents of the fields it writes to 
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=False, n_genes=10)
    
    # get text representation of top DEGs and save to tempdir
    rg = adata.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names  # group names
    lines = []
    for group in groups:
        names = list(rg["names"][group])
        names = adata.var["gene_symbols"][names]
        names = ", ".join(names)
        lines.append(f"DEGs for {group}: {names}")
        if show:
            print(f"DEGs for {group}: {names}")
    with open(os.path.join(output_data_dir, f"DEGs_grouped_by_{groupby}_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.txt"), "w") as f:
        f.write("\n".join(lines))

    # also plot as matrixplot and save to tempdir (CURRENTLY DOES NOT WORK, perhaps too many marker genes?)
    """if marker_file_path != None:
        vprint("Plotting DEG matrixplot...")
        sc.pl.matrixplot(adata, marker_genes_dict, groupby=groupby, dendrogram=False, cmap="viridis", standard_scale="var", colorbar_title="scaled \nexpression", show=False)
        plt.savefig(os.path.join(output_data_dir, f"DEGs_grouped_by_{groupby}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        if show: 
            plt.show() # closes the plot, so run after savefig
            plt.close()"""


def plot_projection_and_DEGs(adata, layer):
    """
    layer should be something like "X" or "X_cnv" or "X_scvi_corrected"

    plots UMAP or PCA projection of most passed obs annotations (skips celltype if isolated) and leiden.
    also computes degs for cleaned obs annotations list.
    """
    # create temp anndata with layer as X, might not have var, if shape of layer and adata.Y don't match
    vprint(f"converting {layer} to anndata...")
    internal_adata = matrix_to_anndata(adata, layer)

    # normalize temp adata and store result in according field in internal adata 
    if not hf.is_normalized(internal_adata) and layer != "X_cnv": # X_cnv is log normalized, so this check does not work (special case)
        vprint("normalizing data...")
        sc.pp.normalize_total(internal_adata, target_sum=1e4)
    else:
        vprint(f"adata is already normalized, skipping to next step...")

    # pca on temp and store result in according field in internal adata
    vprint("computing PCA embeddings...")
    sc.pp.pca(internal_adata, svd_solver="arpack")

    vprint("finding neighbors...")
    sc.pp.neighbors(internal_adata, use_rep="X_pca", n_neighbors=25, metric="euclidean") # these are the default params

    vprint("computing leiden clusters...")
    sc.tl.leiden(internal_adata, resolution=0.1) # uses neighbor graph
    
    # internal list since lists are mutable and we do not want to change the global list
    colored_by = obs_annotations.copy() # smth like ["cell_type", "summed_cnv", "pseudotime", "cnv_score", ...]
    colored_by.append("leiden")
    
    # quality checks before plotting
    vprint("quality checks before plotting...")
    """if "cell_type" in colored_by and cell_type is not None:
        colored_by.remove("cell_type")
        vprint(f"{cell_type} cells are isolated, removing \"cell_type\" from colored_by list")"""
    for entry in colored_by:
        if entry not in internal_adata.obs.keys():
            colored_by.remove(entry)
            vprint(f"{entry} is not in adata.obs, removing it from colored_by list")

    
    # determine suitable number of columns for figure (rougly 1.5 times as many columns as rows)
    n_cols = np.round(np.sqrt(len(colored_by)*1.5)).astype(int) 

    # make sure annotatins with categories are colored by category and not continuous values
    for entry in colored_by:
        if len(internal_adata.obs[entry].unique()) <= 20 or entry == "cnv_clade": # cnv clade yields continuous spectrum prbly bcs it has numpy floats
            vprint(f"turning {entry} into categorical, unique values: {len(internal_adata.obs[entry].unique())}")
            internal_adata.obs[entry] = pd.Categorical(internal_adata.obs[entry])


    # create a figure with one plot per color, save to temp, show if show is True
    if projection == "UMAP":
        vprint("computing UMAP embedding...")
        sc.tl.umap(internal_adata, 0.2) # uses neighbor graph # default value for mindist is actually 0.5 acoording to the docs, not 0.1 (0.2 yielded best seperation of cancer / non cancer)
        sc.pl.umap(internal_adata, color=colored_by, show=False, ncols=n_cols, legend_loc="on data")

        # highlight root cell
        if root_cell_idx != None:
            umap_coords = internal_adata.obsm['X_umap']
            plt.scatter(umap_coords[root_cell_idx, 0], umap_coords[root_cell_idx, 1],
                color='red', s=100, edgecolor='black', label='highlighted cell')

        plt.savefig(os.path.join(output_data_dir, f"UMAP_of_{layer}_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        
        if show: 
            plt.show()
        plt.close()
    elif projection == "PCA":
        sc.pl.pca(internal_adata, color=colored_by, show=False, ncols=n_cols, legend_loc="on data")

        # highlight root cell
        if root_cell_idx != None:
            pca_coords = internal_adata.obsm['X_pca']
            plt.scatter(pca_coords[root_cell_idx, 0], pca_coords[root_cell_idx, 1],
                color='red', s=100, edgecolor='black', label='highlighted cell')

        plt.savefig(os.path.join(output_data_dir, f"PCA_of_{layer}_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        
        if show: 
            plt.show()
        plt.close()
    else:
        raise ValueError(f"projection {projection} is not supported")

    # for each grouping also run deg analysis
    if "projctions+DEGs" in modules:
        if internal_adata.shape[1] == adata.shape[1]: # proxy for varnames being copied over (can't run DEG if you don't know which genes are present)
            
            vprint("Preparping data for DEG analysis...")
            # turn obs annotations into categorical, if they have a reasonable number of categories
            remove_colored_by = []
            for entry in colored_by:
                if len(internal_adata.obs[entry].unique()) <= 50:
                    vprint(f"turning {entry} into categorical, unique values: {len(internal_adata.obs[entry].unique())}")
                    internal_adata.obs[entry] = pd.Categorical(internal_adata.obs[entry])
                else:
                    vprint(f"{entry} has too many categories ({len(internal_adata.obs[entry].unique())}), cannot compute DEGs for it")
                    remove_colored_by.append(entry)

            colored_by = [x for x in colored_by if x not in remove_colored_by] # adjust list

            for entry in colored_by:
                vprint(f"computing DEGs for {entry}")
                get_DEGs(internal_adata, entry)
        else:
            vprint(f"skipping DEG analysis since {layer} has different number of genes than adata")
    else:
        vprint("skipping DEG analysis")


def pseudotime_vs_CNV(adata, show):

    sc.pl.scatter(adata, x="dpt_pseudotime", y="summed_cnvs", title="summed_cnv  vs pseudotime", show = False)
    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"summed_cnvs_vs_pseudotime_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()

    sc.pl.scatter(adata, x="dpt_pseudotime", y="cnv_score", title="cnv_score vs pseudotime", show = False)
    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"cnv_score_vs_pseudotime_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()


def visualize_tree (tree_file, adata = None, obs_columns = None, ring_spacing: float = 0.025, base_radius = 1, target_circumference = 2*np.pi * 0.9, sort_order: Literal["lowest_width_first", "lowest_depth_first"] = "lowest_width_first", show = True, save = False, debug = False):
    """
    Make a circular phylogenetic tree plot. Adata and obs columns used for colored annotation rings around tree.

    Parameters
    ----------

    tree_file : str
        path to tree file (nwk format)
    adata : anndata.AnnData, optional
        anndata object to use (only needed if obs_columns is not empty), by default None
    obs_columns : list, optional
        list of obs columns to use for colored annotation rings, by default None
    ring_spacing : float, optional
        defines distance between annotation rings, by default 0.025 (tight)
    base_radius : float, optional
        base radius of tree plot, by default 1
    target_circumference : float, optional
        target circumference of tree plot, by default 2*np.pi * 0.9 (meaning it will not be a full circle)
    sort_order : Literal["lowest_width_first", "lowest_depth_first"], optional
        defines sort order of leaves, by default "lowest_width_first"
        Currently supports:
            - "lowest_width_first": Clades with the least number of terminal nodes are visited first ("narrowest" clades first).
            - "lowest_depth_first": Clades with the shortest path to a preterminal node (a node that only has terminal nodes as children) are visited first ("shallow" clades first).
            (note that the "lowest_depth_first" option may also choose to visit a generally deep clade that happens to have a single shallow child clade first)
    show : bool, optional
        whether to show the plot, by default True
    save : bool, optional
        whether to save the plot, by default False
    debug : bool, optional
        whether to print debug messages and save debug info to json files in root directory, by default False

    Returns
    -------
    None
    """

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
                    old_name = node.name # debug
                    node.name = f"{str(parents[node.name].name)}_{properly_named_siblings + 1}" # name the node as the nth child of its parent (n = properly_named_siblings + 1, so if a parent has 5 childern, then it will be parent_1 for the first one that is checked, then parent_2, and so on up to root_1_1_1_2_5_1_7_4_2_2...)
                    name_mapping[old_name] = node.name # add new name to mapping

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
        if obs_columns != None:
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
    if type(adata) == sc.AnnData:
        obs = adata.obs.copy()
    else:
        obs = None

    # ASSIGN PRELIMINARY NAMES (needed for functions that work on the tree)
    vprint("Assigning preliminary names to tree nodes...")
    get_preliminary_names(tree)
    # DEBUG: get mapping of preliminary names to actual names
    name_mapping = {node.name: node.name for node in tree.get_nonterminals()} # dict with just preliminary names (values changed in get radial coords)


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
    if obs_columns != None:
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
        for clade in (list(tree.find_clades(order="level"))):
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
        vprint("Showing figure...")
        fig.show()
    if save: 
        print("Oh, you want save to save the figure? That's great! But you will have to do it manually, bcs plotly is a bit weird.")
        print(f"You may want to save to this path: {(os.path.join(output_data_dir, f"{os.path.basename(tree_file).removesuffix('.nwk')}.png"))}")

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

        # DEBUGGING: write name mapping to json
        print("DEBUG: Writing name mapping to json...")
        with open("name_mapping.json", "w") as f:
            json.dump(name_mapping, f, indent=3)
    
    return None


def main():

    # read data into anndata
    print("reading data...")
    print(input_data_file)
    print(type(input_data_file))
    if type(input_data_file) == str and ".h5ad" in input_data_file:
        adata = sc.read_h5ad(input_data_file)
        if verbose:
            vprint("adata summary:")
            vprint(adata) # show summary about adata structure
    else:
        adata = None
        vprint("input_data_file is not a path to a .h5ad file, setting adata to None.")
    
    # limit anndata to cells that fullfill all criteria
    print("limiting cells...")
    adata = limit_cells(adata, selection_criteria).copy()# which column to use and what entry in that column to limit cells to
    print("writing to json (DEBUG)")
    with open (os.path.join(output_data_dir, "selection_criteria.json"), "w") as f: # dump used selection criteria for the run to json
        json.dump(selection_criteria, f, indent=3)


    # UMAP / PCA plots + DEGs for each layer
    if "projections" in modules or "projections+DEGs" in modules:
        for layer in layers: # layerlist like [adata.X, adata.obsm["X_cnv"], ...]
            if layer in adata.layers.keys() or layer in adata.obsm.keys() or layer == "X":
                print(f"plotting projections for {layer}...")
                plot_projection_and_DEGs(adata, layer=layer)
            else:
                raise ValueError(f"layer {layer} not found in adata.layers or adata.obsm or adata.X")
    else:
        vprint("skipping projections module...")

    # scatterplot of pseudotime vs summed_cnvs and cnv_score
    if "pseudotime_vs_cnv" in modules:
        print("plotting pseudotime_vs_cnv...")
        pseudotime_vs_CNV(adata, show=show)
    else:
        vprint("skipping pseudotime_vs_cnv module...")

    # phylogenetic tree
    if "phylogenetic_tree" in modules:
        print("plotting phylogenetic tree...")
        visualize_tree(tree_file, adata, obs_columns=obs_annotations, base_radius=1, target_circumference=target_circumference, sort_order=sort_order, show=show, save=save, debug=False)
    else:
        vprint("skipping phylogenetic_tree module...")

    # OTHER MODULE IDEAS
    # force directed graph for pseudotime
    # DEG heatmap
    # phylogenetic tree for normal / transtion / tumor cells as per cnvscore
    # gene expression over pseudotime graph

    print(f"Output: {output_data_dir}")



if __name__ == "__main__":
    
    # import cmd args
    input_data_file, output_data_dir, modules, selection_criteria, obs_annotations, layers, projection, marker_file_path, root_cell_idx, tree_file, target_circumference, sort_order, show, verbose, save = hf.import_cmd_args(15)
    vprint = hf.make_vprint(verbose)

    if marker_file_path != None:
        with open(marker_file_path, "r") as f:
            marker_genes_dict = json.load(f)
    elif marker_file_path == None and "projections+DEGs" in modules:
        print("No marker file provided, skipping DEG matrix plots...")

    main()

    

