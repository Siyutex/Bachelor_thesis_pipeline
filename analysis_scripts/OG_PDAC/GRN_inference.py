# this script will make a boolean GRN from scRNAseq data (shouled be preprocessed, aggregatedand batch corrected)

import scanpy as sc
import pandas as pd
import numpy as np
import helper_functions as hf
import os
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster


if __name__ == "__main__":
    # consts
    PREFIX = "batch_corrected_HVG_"
    SUFFIX = ".h5ad"


    # get cmd args
    input_data_file, output_data_dir, verbose = hf.import_cmd_args(3)
    vprint = hf.make_vprint(verbose)

    # import data from batch corrected h5ad file
    vprint(f"Importing data from {input_data_file}")
    adata = sc.read_h5ad(input_data_file)

    # check if obs and var names are unique
    if len(adata.obs_names.unique()) != adata.n_obs:
        raise ValueError("Obs names are not unique, obs names should be made unique in aggregate_batches during batch correction.")
    if len(adata.var_names.unique()) != adata.n_vars:
        raise ValueError("Var names are not unique, var names should be unique be default (reference genome should not contain duplicates).")

    # check if obs and var names are unique
    if len(adata.obs_names.unique()) != adata.n_obs:
        raise ValueError("Obs names are not unique, obs names should be made unique in aggregate_batches during batch correction.")
    if len(adata.var_names.unique()) != adata.n_vars:
        raise ValueError("Var names are not unique, var names should be unique be default (reference genome should not contain duplicates).")

    # isolate ductal cells
    vprint("Isolating ductal cells...")
    ductal_cells = adata.obs["cell_type"] == "ductal_cell"
    adata = adata[ductal_cells, :]

    # assign dataframe with var names as column names
    adata_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    print("Dataframe shape: ", adata_df.shape)

    # assign client
    client = Client(LocalCluster())

    # run GRN inference (^2 compute time, 2324 genes take 2:40 minutes, cells do not seem to affect runtime)
    vprint("Running GRN inference...")
    grn = grnboost2(adata_df, verbose=verbose, client_or_address=client, tf_names="all")

    # write GRN to file
    vprint("Writing GRN to file...")
    output_file_path = os.path.join(output_data_dir, f"GRN_{os.path.basename(input_data_file).removesuffix(SUFFIX).removeprefix(PREFIX)}.csv")
    grn.to_csv(output_file_path)

    # send output to executor
    print(f"Output: {output_file_path}")
