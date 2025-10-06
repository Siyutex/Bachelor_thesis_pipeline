# takes edge set from edge rule inference 
# (pseudotime binning, several runs)
# use scikit learn decision tree classifier on each gene + putative regulators to find relevant rules
# turn decision tree into boolean function
# espresso minimization
# export rule set as BNET expressions in JSON, then use for pyboolnet simulation / analysis

from sklearn.tree import DecisionTreeClassifier, export_text
from pyeda.inter import expr, espresso_exprs
import os 
import numpy as np
from scipy.sparse import csc_matrix
import json
import scanpy as sc
import re
import helper_functions as hf

# import cmd args
expression_matrix_file, output_data_dir, verbose, edge_set_file = hf.import_cmd_args(4)
vprint = hf.make_vprint(verbose)

# consts (needed output file naming)
PREFIX = "batch_corrected_HVG_"
SUFFIX = ".h5ad"

def count_bars(line: str):
    """
    Count the number of '|' in a line to determine the depth of that line 
    in a decision tree string.
    """
    return line.count('|')


def parse_condition(line: str) -> str:
    """
    Extract variable and condition from a line in a decision tree string.
    
    Parameters:
        line (str): A line in a decision tree string.
    """

    match = re.search(r"\|--- ([A-Za-z0-9_\-\.]+)\s+(<=|>=|<|>)\s+([0-9.]+)", line) # full alphabet (upper and lower case), all digits, "-" and "." for gene name, then whitespace, then operator, then whitespace, then value
    if not match:
        raise ValueError(f"Could not find a condition in line: {line} in decision tree text")
    var, op, val = match.groups() # assign variables to groups of the match (should be 3 groups, variable, operator, value)

    if op in ('>', '>='):
        vprint(f"Found condition: {var} {op} {val} for line {line}, returning {var} (activator)")
        return var  # greater means the variable is must be 1 for the condition (assuming 0 or 1 are the only values)
    else:
        vprint(f"Found condition: {var} {op} {val} for line {line}, returning ~{var} (inhibitor)")
        return f'~{var}'  # less means the variable is must be 0


def extract_path_rule(lines: list[str], idx: int) -> list[str]:
    """
    Trace back from a line (denoted by idx) containing a leaf node in a decision tree string
    to build a simplified decision tree path, which describes the conditions necessary to
    reach that leaf in the decision tree.

    Should also work to get a path leading to a non-leaf node.

    Parameters:
        lines (list): A full list of lines of a decision tree string.
        idx (int): The index of the line containing the leaf.

    Returns:
        list: A list of conditions that must be met to reach the leaf.
        eg: [A, ~B, C]
    """
    vprint(f"extracting_rule_path for line {lines[idx]}")

    path = []
    current_indent = count_bars(lines[idx])
    
    # loop backwards through the lines until we reach the root (-1 stop needed to actually reach the root at index 0)
    for i in range(idx - 1, -1, -1):
        next_line = lines[i]
        next_indent = count_bars(next_line)

        vprint(f"next_line: {next_line}, next_indent: {next_indent}, current_indent: {current_indent}")

        if next_indent == current_indent - 1: # one level shallower in the decision tree
            condition = parse_condition(next_line) # get the condition of the line (eg A, ~A, ...)
            vprint(f"Condition in next_line: {condition}")
            if condition:
                path.insert(0, condition) # add the condition to the path, insert at the beginning (since we loop backwards)
            current_indent = next_indent
        elif next_indent > current_indent - 1: # skip over deeper or equal levels (they are not part of the path to the leaf node)
            continue
    vprint(f"Extracted path {path}")
    return path


def extract_paths(tree_text: str, target_leaf: str) -> list[str]:
    """
    Find all paths leading to leaf nodes equal to target_leaf and extract 
    rules corresponding to those paths.

    Should also work to get all paths leading to any nodes with a particular value.

    Parameters:
        tree_text (str): A full decision tree string.
        target_leaf (str): The leaf node type to find paths leading to.

    Returns:
        list: A list of lists of conditions that must be met to reach any leave with the 
        value of target_leaf.
        eg: [[A, ~B, C], [A, B, C]]

    """
    vprint(f"Extracting paths leading to {target_leaf}")
    lines = tree_text.strip().split('\n') # get each line in the decision tree
    paths = []
    # get rule for each leaf node that is equal to target_leaf
    for idx, line in enumerate(lines):
        if target_leaf in line:
            path = extract_path_rule(lines, idx)
            if path:
                paths.append(path)
    vprint(f"Paths leading to {target_leaf}:\n{paths}")
    return paths


def simplify_boolean_rule(paths: list[list[str]]) -> str:
    """
    Simplify boolean expression using pyeda/espresso.

    Parameters:
        paths (list): A list of lists of conditions that must be met to reach any leaf.
            (eg: [[A, ~B, C], [A, B, C]], output of extract_paths())

    Returns:
        str: A simplified boolean expression.
    """

    vprint(f"Simplifying rule...")
    # merge paths (AND within a path, OR between paths)
    exprs = [" & ".join(path) for path in paths]
    full_expr = " | ".join(f"({e})" for e in exprs).replace("-","_").replace(".","_DOT_") # encode - and . because pyeda cannot handle them

    parsed_expr = expr(full_expr) # create a pyeda boolean expression that can be used by espresso

    # Use espresso for minimization
    simplified, = espresso_exprs(parsed_expr) # espresso is heuristic, and does not guarantee optimal minimization, but usually gets very close (eg 1 factorization missing). This may affect boolean update time linearly by a factor of maybe 1.2.
    vprint(f"Simplified rule: {simplified}")
    return str(simplified).replace("_","-").replace("_DOT_",".") # decode - and . so they can be used again with the reference genome



def pyeda_to_bnet(expr: str, target: str = "TARGET") -> str:
    """
    Convert a pyeda espresso boolean string into BNET format.
    Example: Or(And(~THY1, COL1A2), And(~THY1, ~CTSK))
    becomes: TARGET, (!THY1 & COL1A2) | (!THY1 & !CTSK)

    IMPORTANT: This assumes disjunctive normal form (eg Or(And(...), And(...)), nested ORs are not supported)

    Parameters:
        expr (str): A pyeda espresso boolean expression.
        target (str, optional): The target node name. Defaults to "TARGET".

    Returns:
        str: A BNET formatted rule.
    """
    # Replace NOT
    expr = expr.replace("~", "!") # looks for ~ and replaces it with !

    # Convert And(...) → ( ... & ... )
    # re.sub takes a regular expression and a replacement function, and a string
    # it then looks for the regex and replaces any matches with the return vlaue of the replacement function
    # then it returns the modified string
    while "And(" in expr:
        expr = re.sub(r"And\(([^()]*)\)",
        lambda m: " & ".join(p.strip() for p in m.group(1).split(",")),
        expr)

    # Convert Or(...) → ( ... | ... )
    while "Or(" in expr:
        expr = re.sub(r"Or\(([^()]*)\)",
        lambda m: " | ".join("(" + p.strip() +")" for p in m.group(1).split(",")),
        expr)

    vprint(f"BNET converted rule: {expr}")
    return f"{target}, {expr}"


# load edge set from json (list of tuples)
vprint(f"Loading edge set from {edge_set_file}")
with open(edge_set_file, "r") as f:
    edge_sets = json.load(f)

# load expression matrix, isolate ductal cells and binarize
vprint(f"Loading expression matrix from {expression_matrix_file}")
adata = sc.read_h5ad(expression_matrix_file)
vprint("Isolating ductal cells...")
ductal_cells = adata.obs["cell_type"] == "ductal_cell"
adata = adata[ductal_cells, :].copy()
binarized_expression_matrix = (adata.X > np.median(adata.X.toarray(), axis=0)).astype(int) # binarize using median (bigger = True = 1, smaller = False = 0)
adata.X = csc_matrix(binarized_expression_matrix)
print(type(binarized_expression_matrix))


# define decsion tree classifier
dtc = DecisionTreeClassifier(criterion="gini",max_depth=5, min_samples_leaf=50) # tune these parameters based on output


# find boolean rules for each edge set
bnet_rules = []
for target, edge_set in edge_sets.items():
    # subset X to training data for current edge set
    vprint(f"Subsetting training data for target {target} and edge set {edge_set}")
    regulator_subset = adata[:, edge_set].copy()
    target_subset = adata[:,adata.var_names == target].copy()

    # train the decision tree classifier
    vprint("Training decision tree classifier...")
    dtc.fit(regulator_subset.X.toarray(), target_subset.X.toarray().flatten())

    # export decision tree as text
    rules = export_text(dtc, feature_names=list(regulator_subset.var_names)).strip()

    # WIP; even if inactive, might be regulator of other, so stil save rule

    if "class: 1" not in rules:
        vprint(f"No no activation function for {target}, all paths lead to deactivation")
        vprint(f"creating constant rule for {target}")
        bnet_rules.append(f"{target}, 0")
    elif "class: 0" not in rules:
        vprint(f"No activation function for {target}, all paths lead to activation")
        vprint(f"creating constant rule for {target}")
        bnet_rules.append(f"{target}, 1")
    elif "class: 1" in rules:
        # Extract all paths leading to "class:1" (means that that path leads to activation)
        paths = extract_paths(rules, "class: 1")

        # Step 2: Simplify Boolean logic with espresso
        simplified_rule = simplify_boolean_rule(paths)

        # Step 3: Convert to BNET format
        bnet_rule = pyeda_to_bnet(simplified_rule, target=target)
        bnet_rules.append(bnet_rule)

# export all rules to json
with open(os.path.join(output_data_dir,f"Boolean_rules_{os.path.basename(expression_matrix_file).removeprefix(PREFIX).removesuffix(SUFFIX)}.json"), "w") as f:
    json.dump(bnet_rules, f)

print("Output: " + os.path.join(output_data_dir, f"Boolean_rules_{os.path.basename(expression_matrix_file).removeprefix(PREFIX).removesuffix(SUFFIX)}.json"))

    











