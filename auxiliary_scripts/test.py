import pandas as pd
import numpy as np
from arboreto.algo import grnboost2
from distributed import Client
from dask import delayed

def main():
    client = Client(n_workers=2, threads_per_worker=4, memory_limit="4GB")
    print(client)

    data = pd.DataFrame({
        f"gene{i}": [i + np.random.normal(100, 50) for _ in range(1000)]
        for i in range(20)  # test with 20 genes first
    }, index=[f"cell{j}" for j in range(1000)])

    @delayed
    def add(a, b):
        return a + b

    c = add(2,3)
    print(c.compute())

    network = grnboost2(data, client_or_address=client, tf_names="all", verbose=True)
    print(network.head())

if __name__ == "__main__":
    main()