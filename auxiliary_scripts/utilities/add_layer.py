import scanpy as sc

FILE_LOCATION = r"/home/julian/Bachelor_thesis_pipeline/Data/output_storage/isolated/isolated_transition_clades_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected.h5ad"

def add_log1p_layer(adata):
    adata_log = sc.pp.log1p(adata, copy=True)
    adata.layers["log1p"] = adata_log.X


if __name__ == "__main__":
    print("reading adata...")
    adata = sc.read_h5ad(FILE_LOCATION)

    print("Adding log1p layer...")
    add_log1p_layer(adata)

    print("Saving results...")
    adata.write(FILE_LOCATION, compression="gzip")