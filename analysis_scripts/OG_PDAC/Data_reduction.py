import scanpy as sc
import numpy as np
import helper_functions as hf
import os


def reduce_adata(adata, layers_to_remove):
    for element in layers_to_remove:
        if element in adata.obsm.keys():
            del adata.obsm[element]
        elif element in adata.layers.keys():
            del adata.layers[element]
        elif element == "X":
            del adata.X
        else:
            vprint(f"Element {element} not found in adata.obsm or adata.layers. Cannot remove.")

def main():

    # load adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)
    vprint(f"Adata summary:\n{adata}")

    # reduce obsm and layers
    if layers_to_remove:
        print("reducing layers and obsm matrices...")
        reduce_adata(adata,layers_to_remove)
    vprint(f"Adata summary after reduction:\n{adata}")

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))


if __name__ == "__main__":

    input_data_file, output_data_dir, layers_to_remove, verbose = hf.import_cmd_args(4)
    vprint = hf.make_vprint(verbose)

    main()