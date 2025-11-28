import scanpy as sc
import helper_functions as hf


def get_mito_percentage_by_cnvscore(adata, layer, cnv_score):
    # print average mito percentage per cell with cnv score > cnv_score
    
    # isolate cells with cnv score > cnv_score
    internal_adata = hf.matrix_to_anndata(adata, layer)
    adata_sub = internal_adata[internal_adata.obs["cnv_score"] > cnv_score]

    # get average mito percentage per cell
    mito_percentage = adata_sub.obs["pct_counts_mito"].mean()

    print(f"Amount of cells with cnv_score > {cnv_score}: {adata_sub.shape[0]}")
    print(f"Average mitochondrial gene expression in cells with cnv_score > {cnv_score}: {mito_percentage:.3f}")


def get_global_highest_mito_percentage(adata):
    # find cell with highest pct_counts_mito

    mito_percentage = adata.obs["pct_counts_mito"].max()

    print(f"Cell with highest mitochondrial percentage: {adata.obs[adata.obs['pct_counts_mito'] == mito_percentage].index[0]}")
    print(f"Mitochondrial percentage of this cell: {mito_percentage:.3f}")

if __name__ == "__main__":

    FILE_LOCATION = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/preprocessed/preprocessed_PDAC_cancerous_0.h5ad"

    print("Reading data...")
    adata = sc.read_h5ad(FILE_LOCATION)

    # print("Getting mitochondrial percentage...")
    # get_mito_percentage_by_cnvscore(adata, "X_scANVI_corrected", 0.01)

    print("Finding cell with highest mitochondrial expression")
    get_global_highest_mito_percentage(adata)