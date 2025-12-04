import scanpy as sc
import os

def remove_annoation(adata, annotation_key):
    if annotation_key in adata.obs.keys():
        del adata.obs[annotation_key]
    if annotation_key in adata.var.keys():
        del adata.var[annotation_key]
    if annotation_key in adata.obsm.keys():
        del adata.obsm[annotation_key]
    if annotation_key in adata.layers.keys():
        del adata.layers[annotation_key]
    if annotation_key == "X":
        del adata.X

    return adata


if __name__ == "__main__":
    input_file = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN3.6/tree/transition_clades_PDAC_ductal_cell.h5ad"
    output_dir = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN3.6/reduced"
    adata = sc.read_h5ad(input_file)
    for annotation in ["cnv_clade", "cancer_state_inferred_tree"]:
        adata = remove_annoation(adata, annotation)

    print(f"Adata after removal of annotations:\n{adata}")

    adata.write(os.path.join(output_dir, os.path.basename(input_file)), compression="gzip")