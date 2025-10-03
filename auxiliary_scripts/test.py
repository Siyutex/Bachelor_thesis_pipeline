import numpy as np

subsample_idxs = np.random.randint(0, 1000, 1000)
print(subsample_idxs)

print(len(list(set(subsample_idxs))))
