import os
import scanpy as sc
import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import resample
from sklearn.preprocessing import RobustScaler
import helper_functions as hf

def check_lib_size(adata, layer_to_check):
    if layer_to_check in adata.layers.keys():
        print(f"Library size of {layer_to_check} is {adata.layers[layer_to_check].sum(axis=1)}")
    elif layer_to_check in adata.obsm.keys():
        print(f"Library size of {layer_to_check} is {adata.obsm[layer_to_check].sum(axis=1)}")

    if adata.X is not None:
        print(f"Head of adata.X: \n{adata.X[:5]}")
    for layer in adata.layers.keys():
        print(f"Layer {layer} summary: \n{adata.layers[layer]}")
    for layer in adata.obsm.keys():
        print(f"Obsm {layer} summary: \n{adata.obsm[layer]}")


def show_annotation(adata, extract):
    # check a particular obs / var annotation (type, first few entries)
    if extract in adata.obs.keys():
        print(f"Type of entries in {extract}: {type(adata.obs[extract][0])}")
        print(f"First few entries of {extract}: \n{adata.obs[extract][:5]}")
    elif extract in adata.var.keys():
        print(f"Type of entries in {extract}: {type(adata.var[extract][0])}")
        print(f"First few entries of {extract}: \n{adata.var[extract][:5]}")


def check_n_obs(adata):
    # check how many observations there are
    print(f"Number of observations in X / layers: {adata.n_obs}")
    for key in adata.obsm.keys():
        print(f"Number of observations in {key}: {adata.obsm[key].shape[0]}")

def check_obs_percentage(adata, obs_dict):
    """
    Check what percentage of cells has obs_column == value
    where obs column is they key in the dict and value is the union of the values
    """

    for obs_column, value in obs_dict.items():
        mask = adata.obs[obs_column].isin(value)
        print(f"Percentage of cells with {obs_column} == {value}: {mask.sum() / adata.n_obs * 100}")


def calculate_hopkins_stable(adata, n_comps=50, m=1000, iterations=10):
    """
    Calculates a robust Hopkins Statistic for an AnnData object.
    Uses cells from IQR to avoid outliers affecting bounding box sizes for random points.
    
    H ~ 1.0 : Highly Clustered (Low Noise/Clear States)
    H ~ 0.5 : Random/Uniform (High Noise/Blurred Manifold)
    
    Parameters:
    -----------
    adata : sc.AnnData
        The annotated data matrix.
    n_comps : int
        Number of principal components to use for dimensionality reduction.
    m : int
        Number of cells to sample (keeps density constant across datasets).
    iterations : int
        Number of resamplings to perform for a stable estimate.
    """

    adata = sc.read_h5ad(path)
    adata = hf.matrix_to_anndata(adata, matrix_key="X_scANVI_corrected").copy()
    sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack")
    
    # Extract the high-dimensional embedding (e.g., 30D PCA)
    X = adata.obsm["X_pca"]
    n_cells, n_dims = X.shape
    
    # Robust Scaling: Center and scale based on percentiles.
    # This prevents extreme outliers from stretching the space.
    # It also makes it so that each dimension of the data contributes evenly to the distance (because it scales each feature seperately).
    # this means the 1 or 2 extremely clustered dimensions (eg PC 1 and 2) don'T dominate the statistic and make a noisy dataset seem clustered (it makes the metric stricter).
    # But if there are 1 or 2 extremely clustered dimensions, the real points will still take much less volume in the bounding box and this will increase H.
    X_scaled = RobustScaler().fit_transform(X)
    
    hopkins_scores = []
    
    for i in range(iterations):
        # Fit the Nearest Neighbors model once per sample (to make sure density is equal)
        # We use n_neighbors=2 because the 1st neighbor to a real cell is itself
        neigh = NearestNeighbors(n_neighbors=2).fit(X_scaled)

        # Subsample 'm' real cells (ensures density consistency)
        # eg if you were to sample the same dataset / similar datasets twice with different densities, then it would be more likely that there are some random real points spread throughout -> lower H.
        X_sample = resample(X_scaled, n_samples=m, replace=False)
        
        # Create the "Bounding Box" using 5th and 95th percentiles  (again to avoid outliers)
        # to ensure random points live where the biological signal is.
        mins = np.percentile(X_scaled, 5, axis=0)
        maxs = np.percentile(X_scaled, 95, axis=0)
        
        # Generate random points in the same n-dimensional space
        # np.random.uniform uses the vectors 'mins' and 'maxs' to define the box.
        random_points = np.random.uniform(mins, maxs, (m, n_dims))
        
        # Distances from random points to nearest real data (u)
        # compares dataset used to initialize the object to the dataset passed
        # so in this case it returns the distance, index of the first nearest neighbour in the real dataset for each point in the random points
        u_dist, _ = neigh.kneighbors(random_points, n_neighbors=1)
        u_sum = np.sum(u_dist)
        
        # Distances from real sampled points to their nearest real neighbor (w)
        w_dist, _ = neigh.kneighbors(X_sample, n_neighbors=2)
        w_sum = np.sum(w_dist[:, 1]) # Index 1 is the actual neighbor
        
        # Calculate H for this iteration
        H = u_sum / (u_sum + w_sum)
        hopkins_scores.append(H)
    
    # Return the median across all iterations for stability
    print(f"Median Hopkins statistic: {np.median(hopkins_scores)}")

if __name__ == "__main__":

    layer_to_check = None # check library sizes for this layer
    dir = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/subsampled/PDAC/original_subsamples"

    obs_dict = {
        "cancer_state_inferred_tree": ["transitional"],
    }

    for file in os.listdir(dir):
        path = os.path.join(dir, file)    
        check_n_obs(sc.read_h5ad(path))

        



