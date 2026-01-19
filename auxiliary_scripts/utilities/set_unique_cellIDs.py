# takes a mtx files as input, loads into anndata, sets unique cell IDs, saves to mtx file
# USAGE: 
# - fill in output_dir
# - add directory paths with files to dir_list
# - set var_names 
import scanpy as sc
import os


if __name__ == "__main__":
    # paramteters
    OUTPUT_DIR = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/OG_data/Shin_et_al./manual_cell_IDs_shin"
    DIR_LIST = [
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/OG_data/Shin_et_al./shin_cancerous",
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/OG_data/Shin_et_al./shin_non_cancerous",
    ]
    VAR_NAMES = "gene_ids" # gene_ids or gene_symbols



    # script
    file_list = [os.path.join(dir, file) for dir in DIR_LIST for file in os.listdir(dir)]

    start_id = 0
    # load file, make unique cell IDs, save file, repeat
    for file_path in file_list:
        # load adata and save number of cells
        print(f"reading {file_path}")
        adata = sc.read_10x_mtx(file_path, var_names=VAR_NAMES) # manually set gene_ids for ENSG IDs or gene_symbols
        n_cells = adata.n_obs

        # set unique cell IDs
        unique_obs_names = [f"cell_{i+start_id}" for i in range(len(adata.obs_names))]
        adata.obs_names = unique_obs_names

        # update start_id (so names stay unique across files)
        start_id += n_cells

        # save adata
        final_output_path = os.path.join(OUTPUT_DIR, os.path.basename(file_path) + ".h5ad")
        print(f"writing {final_output_path}")
        adata.write(final_output_path, compression="gzip")