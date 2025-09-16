import numpy as np

def make_vprint(verbose):
    """ 
    Return a function that prints if verbose, the argument of this function, is true.
    This is supposed to be called once at the beginning of a function definition.

    Parameters
    ----------
    verbose : bool
        Whether to print the output or not

    Returns
    -------
    function
        A lambda function that prints if verbose is true

    Example
    -------
    >>> vprint = make_vprint(verbose)
    >>> vprint("Hello world")
    """
    return lambda *a, **k: print(*a, **k) if verbose else None 
    # lambda notation: lambda argument1, arugment2: return value (so here we return a function that prints if verbose is true)
    # *a is a tuple of positional arguments, **k is a dictionary of keyword arguments
    # eg in function(a,b,c, d=1, e=2) *a = (a,b,c) and **k = {"d":1, "e":2}

# check if data is alrdy normalized
def is_normalized(adata, layer: str = None, verbose: bool = False):

    """
    Check if adata is already normalized. Checks for the given layer, which should be a matrix of shape adata.X.shape.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data object.
    layer : str, optional
        Which layer to check for normalization. If None, defaults to adata.X.
    verbose : bool, optional
        Whether to print additional information. Defaults to False.

    Returns
    -------
    bool
        Whether the data is normalized.

    Notes
    -----
    Checks if the sum of counts per cell is the same for all cells in the layer.
    If verbose, prints the number of barcodes, median UMIs per barcode, and maximum UMIs per barcode.
    """
    
    vprint = make_vprint(verbose) # function to print if verbose

    # sum of counts per cell for first 100 cells
    if layer == None:
        counts_per_cell = np.array(adata.X[:100].sum(axis=1)).flatten()
    else:
        counts_per_cell = np.array(adata.layers[layer][:100].sum(axis=1)).flatten()

    if verbose == True:
        vprint(f"   Number of barcodes: {adata.n_obs}")
        vprint(f"   Median UMIs per barcode: {np.median(counts_per_cell):.1f}")
        vprint(f"   Max UMIs per barcode: {np.max(counts_per_cell):.0f}")

    # check if all counts are equal
    if np.allclose(counts_per_cell, counts_per_cell[0]): # verify that all entries are equal (integer precision, not useful for comparing floats)
        if verbose == True:
            vprint("adata is normalized")
        return True
    else:
        if verbose == True:
            vprint("adata is NOT normalized")
            vprint(f"Sums of counts per cell for first 100 rows of adata.X:\n{counts_per_cell}")
        return False
