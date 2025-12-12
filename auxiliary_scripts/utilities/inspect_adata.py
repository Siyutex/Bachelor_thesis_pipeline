import scanpy as sc


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

if __name__ == "__main__":

    layer_to_check = None # check library sizes for this layer
    file_locations = [
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/scMF/scMF_CNV_inferred_shin.h5ad",
    ]



    obs_dict = {
        "cancer_state": ["cancerous"],
        "cancer_state_inferred": ["cancerous"],
        "cancer_state_inferred_scMF": ["cancerous"]
    }

    for path in file_locations:
        adata = sc.read_h5ad(path)
        check_obs_percentage(adata, obs_dict=obs_dict)
        



