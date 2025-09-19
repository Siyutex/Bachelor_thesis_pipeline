# subprocess of Pipeline_executor.py
# input should be a preprocessed h5ad file (can be compressed)
# uses leiden clustering in PCA space + known markers for each tissue type to annotate cell type
# different treatment for cancer cells, bcs dedifferentiation

import scanpy as sc
import os
import sys
from scipy.sparse import issparse
import re
import json

# import command line arguments from ececutor script
input_data_file = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input file")
output_data_dir = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")

# import input data into anndata object (handles compressed h5ad files as well)
adata = sc.read_h5ad(input_data_file)

#-------------------------------------
# clustering and associated processing
#-------------------------------------
def process_and_cluster(adata):


    # normalization (needed to account for cells with different total expression levels)
    """
    Process and cluster anndata object.

    Parameters
    ----------
    adata : AnnData
        Anndata object to process and cluster.

    Returns
    -------
    adata : AnnData
        Anndata object with processed data and leiden cluster labels stored in adata.obs['leiden' | key_added].
    batch : str
        String indicating the batch of origin
    cancer_state : str
        String indicating the cancer state of the tissue type (e.g. 'cancerous', 'non_cancerous').
    """
    print("Filtering mitochondiral and ribosomal genes before DEG analysis")
    ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))
    mito_genes = adata.var_names.str.startswith('MT-')
    adata = adata[:, ~(ribo_genes | mito_genes)]

    print("normalizing data")
    sc.pp.normalize_total(adata, target_sum=1e4)

    # logarithmization (needed to reduce skew from highly expressed genes)
    print("logarithmizing data")
    sc.pp.log1p(adata)

    # DO NOT SCALE expression, how would you find highly expressed genes per cluster if all genes are centered at 0 mean?

    # dimensionality reduction
    print("computing PCA")
    sc.pp.pca(adata, svd_solver="arpack") # ads the PCA representation to adata.obsm['X_pca' | key_added]

    # generate k nearest neighbor graph (required for leiden clustering)
    print("computing neighbors")
    sc.pp.neighbors(adata, n_neighbors=30) # consider twice the normal amount of neighbors for bigger clusters

    # assign cluster labels
    print("computing leiden clusters")
    sc.tl.leiden(adata) # adds the leiden cluster labels to adata.obs['leiden' | key_added]
    # print number of cluster and number of cells per cluster
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
    print(adata.obs['leiden'].value_counts())

    # get tissue type and cancer state (also contains cancer state)
    batch = os.path.basename(input_data_file).removeprefix("preprocessed_").removesuffix(".h5ad")
    if "non_cancerous" in batch:
        cancer_state = "nc"
    elif "cancerous" in batch:
        cancer_state = "c"

    print(batch)
    print(cancer_state)

    return adata, batch, cancer_state

#-------------------------------------
# cell type annotation
#-------------------------------------


# loop through each cluster
# check which genes rank highest in expression in that cluster
# compare to known markers for current tissue type (use current tissue type as function argument)
# annotate anndata with cell type for each cluster (if leiden = A then cell type = X)
# export as h5ad with gz compression


def annotate_cell_type(adata, batch, cancer_state):
    """
    Annotates the cell type for each cluster in the data.

    Parameters
    ----------
    adata : AnnData
        The input data
    tissue_type : str
        The tissue type

    Returns
    -------
    adata : AnnData
        The annotated data
    """

    # NOTE: cannot just take most highly expressed genes, because each cell has a bunch of actin and rsp and so on
    # need sig diff genes per cluster

    # run differential expression across clusters
    print("computing differentially expressed genes")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=5, use_raw=False)  # limit to 5 genes and don't use raw matrix for speed (would not make sense to process matrix and then not use it), 6 jobs, one per CPU core

    # top DEGs are stored in adata.uns["rank_genes_groups"]["names"][groups / clusters][top DEG list]
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    # get top 5 marker genes for each cluster
    print("getting top 5 marker genes for each cluster")
    top_markers = {}
    for group in groups:
        top_markers[group] = result["names"][group][:5].tolist()
        print(f"Top {len(top_markers[group])} genes in cluster {group}: {top_markers[group]}")#
    
    return top_markers




        

adata, batch, cancer_state = process_and_cluster(adata)
top_markers = annotate_cell_type(adata, batch, cancer_state)

# export marker dictionary as json
with open(os.path.join(output_data_dir, f"{batch}_markers.json"), "w") as f:
    json.dump(top_markers, f)

print("Output: " + os.path.join(output_data_dir, f"{batch}_markers.json"))