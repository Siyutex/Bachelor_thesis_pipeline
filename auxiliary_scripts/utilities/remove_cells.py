import scanpy as sc
import numpy as np


def remove_high_mito(adata, mito_percentage):
    print(f"Amount of cells before high mito removal: {adata.shape[0]}")

    if adata.obs.get("pct_counts_mito") is None:
        if adata.var_names.str.startswith("ENSG").all():
            adata.var["mito"] = adata.var["gene_symbols"].str.startswith("MT-")  # identify mitochondrial genes, assuming they start with "MT-"
        else:
            adata.var["mito"] = adata.var_names.str.startswith("MT-")
        adata.obs["pct_counts_mito"] = adata.X[:, adata.var["mito"].values].sum(axis=1) / adata.X.sum(axis=1)
    
    adata = adata[adata.obs['pct_counts_mito'] < mito_percentage, :]

    print(f"Amount of cells after high mito removal: {adata.shape[0]}")


if __name__ == "__main__":
    input_data_file = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/batch_corrected/batch_corrected_PDAC.h5ad"

    print("reading data")
    adata = sc.read_h5ad(input_data_file)
    print("removing high mito cells")
    remove_high_mito(adata, 0.15)
    print("saving data")
    adata.write(input_data_file.removesuffix(".h5ad")+"_mito_removed.h5ad", compression="gzip")