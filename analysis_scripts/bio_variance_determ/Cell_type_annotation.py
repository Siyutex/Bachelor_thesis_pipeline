# subprocess of Pipeline_executor.py
# input should be a preprocessed h5ad file (can be compressed)
# uses leiden clustering in PCA space + known markers for each tissue type to annotate cell type
# different treatment for cancer cells, bcs dedifferentiation

import scanpy as sc
import os
import sys
import tempfile
import numpy as np
from scipy.sparse import issparse

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
    tissue_type : str
        String indicating the tissue type (e.g. 'colon','glioma' , etc.).
    cancer_state : str
        String indicating the cancer state of the tissue type (e.g. 'cancerous', 'non_cancerous').
    """
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
    sc.pp.neighbors(adata)

    # assign cluster labels
    print("computing leiden clusters")
    sc.tl.leiden(adata) # adds the leiden cluster labels to adata.obs['leiden' | key_added]

    # get tissue type (also contains cancer state)
    tissue_type = os.path.basename(input_data_file).removeprefix("preprocessed_").removesuffix(".h5ad")
    if "non_cancerous" in tissue_type:
        cancer_state = "c"
        tissue_type = tissue_type.removesuffix("_non_cancerous")
    elif "cancerous" in tissue_type:
        cancer_state = "nc"
        tissue_type = tissue_type.removesuffix("_cancerous")

    return adata, tissue_type, cancer_state

#-------------------------------------
# cell type annotation
#-------------------------------------


# loop through each cluster
# check which genes rank highest in expression in that cluster
# compare to known markers for current tissue type (use current tissue type as function argument)
# annotate anndata with cell type for each cluster (if leiden = A then cell type = X)
# export as h5ad with gz compression


def annotate_cell_type(adata, tissue_type, cancer_state):
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

    # define known marker dictionary
    marker_dict = {
        'colon_c': [],
        'colon_nc': [],
        'glioma_c': [],
        'glioma_nc': [],
        'breast_cancer_c': [],
        'breast_cancer_nc': [],
        'CCRC_c': [],
        'CCRC_nc': [],
    }

    # loop through each cluster
    for cluster in adata.obs['leiden'].unique():
        # subset of adata for the cluster
        cluster_cells = adata[adata.obs['leiden'] == cluster].X

        mean_exp = np.array(cluster_cells.X.mean(axis=0)).flatten()
        top_5_idx = mean_exp.argsort()[::-1][:5] # [::-1] reverses the array, [:5] takes the first 5 elements
        top_5_genes = adata.var_names[top_5_idx] # list of top 5 gene names in the cluster (should be enough to identify cell type)


        # now go through entire cellmarker 2 database (xlsx) 
        # for every marker append corresponding cell name to a list
        # the most common cell name is the annotation, if equal amounts then unknown
        # also put this out as debugging info
        # then based on first few datasets curate a manual list of cell types and corresponding markers
        # annote all ccrc datasets
        # scANVI
        # isi ai schisgamber di is squa'im dusch kob labam
        
            
