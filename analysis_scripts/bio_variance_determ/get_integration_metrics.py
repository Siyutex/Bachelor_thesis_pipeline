# This script takes a batch corrected h5ad file as input
# It should have bathces, cancer states (or similar conditions, like cell types) annotated in adata.obs
# It should have a corrected expression matrix in adata.obsm
# These go into the global constants
# The script computes integration metrics (currently kBET only), and DEG overlap between the corrected and raw matrices
# it imports functions from global_scripts. As such, it has to be run from the command line in the root of the repository
# to do this run python -m analysis_scripts.bio_variance_determ.get_integration_metrics

import scanpy as sc
import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import chisquare
from global_scripts import helper_functions as hf

# ---- Inputs ----
DATA_INPUT = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\debugging_tmps\test.h5ad"
adata = sc.read_h5ad(DATA_INPUT)
print(f"Loaded adata from {DATA_INPUT}")

# ---- Constants ----
# adata adjacent constants and keys
BATCH_KEY = "batch" #in adata.obs[BATCH_KEY], per cell
CONDITION_KEY = "cancer_state" # in adata.obs[CONDITION_KEY], per cell
CORRECTED_MATRIX_KEY = "X_scVI_corrected" # in adata.obsm[CORRECTED_MATRIX_KEY], per integrated dataset

# script relevant constants
N_REDUCED_DIMENSIONS = 10 # number of dimensions to reduce to (eg in PCA)
NORMALIZATION_VALUE = 10000 # counts per cell to normalize to


# ---- Functions ----
#kBET implementation
def compute_kBET(adata, embedding: str = None, n_neighbors: int = 15, alpha: float = 0.05, verbose: bool = False):
    vprint = hf.make_vprint(verbose)
    
    ad = adata.copy()  # keep original safe

    vprint("Computing neighbors")
    sc.pp.neighbors(ad, use_rep=embedding, n_neighbors=n_neighbors)

    vprint("Computing global batch proportions")
    batches = ad.obs["batch"].to_numpy()
    labels, counts = np.unique(batches, return_counts=True)
    global_props = dict(zip(labels, counts / counts.sum())) # put global batch proportions into a dict, key = batch name, value = proportion(%)
    vprint(f"Global batch proportions: \n {global_props}")

    vprint("Finding batch identities of k nearest neighbors")
    dist = ad.obsp["distances"].tocsr() # info on csr matrix: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
    neighbors_batches = [
        batches[dist.indices[dist.indptr[i]+1:dist.indptr[i+1]]]
        for i in range(dist.shape[0])
    ]

    vprint("Running chi2 test for each cell + neighbors")
    p_values = []
    for nb in neighbors_batches:
        observed = Counter(nb)
        expected = np.array([global_props[k] * len(nb) for k in labels])
        observed_arr = np.array([observed.get(k, 0) for k in labels])
        _, pval = chisquare(f_obs=observed_arr, f_exp=expected) # chi square statistic is sum of (observed - expected)^2 / expected for each batch. larger if observed and expected counts differ (eg I expect 3 from batch 1 but get 7). Then compare to multinomial distribution to get p value
        p_values.append(pval)

    rejection_rate = np.mean(np.array(p_values) < alpha)
    return rejection_rate


def compute_DEG_overlap(adata, verbose: bool = False):
    vprint = hf.make_vprint(verbose)

    # Run DE on raw log-normalized counts
    vprint("computing DEGs for raw counts")
    sc.tl.rank_genes_groups(adata, groupby=CONDITION_KEY, method="wilcoxon", use_raw=False, n_genes=200)
    genes_pre = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).iloc[:,0].values[:200]

    # Run DE on scVI-corrected expression (X_scVI_corrected)
    vprint("computing DEGs for corrected expression")
    sc.tl.rank_genes_groups(adata, groupby=CONDITION_KEY, method="wilcoxon", layer=CORRECTED_MATRIX_KEY, n_genes=200)
    genes_post = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).iloc[:,0].values[:200]

    overlap = len(set(genes_pre) & set(genes_post)) / 200

    return overlap


def compute_metrics(adata, embed, verbose=False):
    vprint = hf.make_vprint(verbose)

    metrics = {}
    # kBET (are the batches mixed well locally / within clusters, 0 = perfect mixing, 1 = no mixing)
    vprint(f"Computing kBET for {embed}")
    metrics["kBET"] = compute_kBET(adata, embedding=embed, verbose=verbose)

    return metrics









# ---- Validation ----
# Check that all necessary keys are present
assert BATCH_KEY in adata.obs
assert CONDITION_KEY in adata.obs
assert CORRECTED_MATRIX_KEY in adata.obsm


# ---- Preprocessing ----
# Add a PCA embedding from raw counts for comparison, check normalization
if not hf.is_normalized(adata):
    sc.pp.normalize_total(adata, target_sum=NORMALIZATION_VALUE)
sc.pp.pca(adata, n_comps=N_REDUCED_DIMENSIONS) # 10 components so identical to scVI which has 10 latent dimensions

# Also add a PCA embedding for the corrected data
if CORRECTED_MATRIX_KEY in adata.obsm:
    adata.layers[CORRECTED_MATRIX_KEY] = adata.obsm[CORRECTED_MATRIX_KEY]
    if not hf.is_normalized(adata, layer=CORRECTED_MATRIX_KEY):
        sc.pp.normalize_total(adata, layer=CORRECTED_MATRIX_KEY, target_sum=NORMALIZATION_VALUE)
    sc.pp.pca(adata, layer=CORRECTED_MATRIX_KEY, n_comps=N_REDUCED_DIMENSIONS, key_added=CORRECTED_MATRIX_KEY + "_pca")
else:
    raise KeyError("No corrected expression matrix found for the chosen CORRECTED_MATRIX_KEY. Please provide one in adata.obsm[CORRECTED_MATRIX_KEY]")


# ---- Metric computation ----
results = {}
for embed in ["X_pca", CORRECTED_MATRIX_KEY + "_pca"]:
    results[embed] = compute_metrics(adata, embed, verbose=True)

df = pd.DataFrame(results).T

# DEG overlap (are the distinguishing genes per condition (in CONDITION_KEY) the same before and after integration, 0 = no overlap, 1 = perfect overlap)
print(f"Computing DEG overlap for {embed}")
DEG_overlap = compute_DEG_overlap(adata, verbose=True)

df["DEG_overlap"] = DEG_overlap
print(df)
