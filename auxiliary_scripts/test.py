from Bio.Phylo.BaseTree import Clade, Tree
import random

def random_binary_tree(n_leaves, name_prefix="T"):
    """Generate a random binary tree with n_leaves terminal nodes."""
    # start with leaf clades
    leaves = [Clade(name=f"{name_prefix}{i+1}") for i in range(n_leaves)]
    
    # keep combining randomly until one root remains
    while len(leaves) > 1:
        # randomly choose two clades to merge
        a, b = random.sample(leaves, 2)
        leaves.remove(a)
        leaves.remove(b)
        # create a new parent clade
        new_clade = Clade()
        new_clade.clades.extend([a, b])
        leaves.append(new_clade)
    
    # the last remaining clade is the root
    return Tree(root=leaves[0])

# Example usage
tree = random_binary_tree(20)
from Bio import Phylo
Phylo.draw_ascii(tree)

# write the tree to a Newick file
Phylo.write(tree, "50TNtree.nwk", "newick")