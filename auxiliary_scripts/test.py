# for a given decision tree like this
"""
Subprocess GRN_rule_inference.py: |--- THY1 <= 0.50

Subprocess GRN_rule_inference.py: |   |--- COL1A2 <= 0.50

Subprocess GRN_rule_inference.py: |   |   |--- CTSK <= 0.50

Subprocess GRN_rule_inference.py: |   |   |   |--- class: 0

Subprocess GRN_rule_inference.py: |   |   |--- CTSK >  0.50

Subprocess GRN_rule_inference.py: |   |   |   |--- class: 0

Subprocess GRN_rule_inference.py: |   |--- COL1A2 >  0.50

Subprocess GRN_rule_inference.py: |   |   |--- CTSK <= 0.50

Subprocess GRN_rule_inference.py: |   |   |   |--- class: 0

Subprocess GRN_rule_inference.py: |   |   |--- CTSK >  0.50

Subprocess GRN_rule_inference.py: |   |   |   |--- class: 0

Subprocess GRN_rule_inference.py: |--- THY1 >  0.50

Subprocess GRN_rule_inference.py: |   |--- class: 1
"""
# in order to extract all rules traverse the char array and find all indexes of "1"
# then for each such index, traverse the char array backwards from that index (first iteration at lowest index, then ascending order) and extract a simplified tree like:
"""
    Subprocess GRN_rule_inference.py: |--- THY1 >  0.50

    Subprocess GRN_rule_inference.py: |   |--- class: 1
    """
    # this can be done by traversing backwards until a an isolated instance of "|" is found (means we are at the root of the tree)
    # if the current lowest count of "|" is n and in a previous line the count is n+1, skip that line
    # then just export that string to a list
# now we have a list of simplified trees
# then turn each tree into a pyboolnet rule (>= or > means, must be 1 (ie A), < or <= mean must be 0 (ie ~A))
# connect all rules within one tree with AND
# then connect all trees with OR
# then run espresso logic minimization
# export simplified rule for each target gene to JSON




import re
import json
from pyeda.inter import expr, espresso_exprs


def count_bars(line):
    """Count the number of '|' in a line to determine the depth."""
    return line.count('|')


def parse_condition(line):
    """Extract variable and condition from a line."""
    print(f"Extracting condition from line: {line}")
    match = re.search(r"\|--- ([A-Za-z0-9_\-]+)\s+(<=|>=|<|>)\s+([0-9.]+)", line)
    print(f"Match: {match.groups()}")
    if not match:
        raise ValueError(f"Could not parse condition from line: {line} in decision tree text")
    var, op, val = match.groups()
    if op in ('>', '>='):
        return var  # A = 1
    else:
        return f'~{var}'  # A = 0


def extract_rule_path(lines, idx):
    """Trace back from a class: 1 line to build a simplified path."""
    path = []
    current_indent = count_bars(lines[idx])
    print(f"extract_rule_path, current_indent: {current_indent} in line: {lines[idx]} (idx: {idx})")
    
    for i in range(idx - 1, -1, -1):
        line = lines[i]
        indent = count_bars(line)

        if indent == current_indent - 1:
            condition = parse_condition(line)
            if condition:
                path.insert(0, condition)
                print(f"Added condition: {condition}")
            current_indent = indent
        elif indent > current_indent - 1:
            # skip over deeper levels
            continue
    print(f"Extracted path {path}")
    return path


def extract_all_class1_paths(tree_text):
    print(f"extract_all_class1_paths, tree_text: {tree_text}")
    """Find all class: 1 lines and extract rule paths for each."""
    lines = tree_text.strip().split('\n')
    print(f"extract_all_class1_paths, lines: {lines}")
    paths = []
    for idx, line in enumerate(lines):
        if "class: 1" in line:
            path = extract_rule_path(lines, idx)
            if path:
                paths.append(path)
    return paths


def simplify_boolean_rule(paths):
    """Simplify boolean expression using pyeda/espresso."""
    print(f"Simplifiy_boolean_rule, paths: {paths}")
    # Create OR of ANDs
    exprs = [" & ".join(path) for path in paths]
    print(f"Simplifiy_boolean_rule, exprs: {exprs}")
    full_expr = " | ".join(f"({e})" for e in exprs).replace("-","_")
    print(f"Simplifiy_boolean_rule, full_expr: {full_expr}")
    parsed_expr = expr(full_expr)
    print(f"Simplifiy_boolean_rule, parsed_expr: {parsed_expr}")

    # Use espresso for minimization
    simplified, = espresso_exprs(parsed_expr) # espresso is heuristic, and does not guarantee optimal minimization, but usually gets very close (eg 1 factorization missing). This may affect boolean update time linearly by a factor of maybe 1.2.
    return str(simplified).replace("_","-")



def pyeda_to_bnet(expr: str, target: str = "MYNODE") -> str:
    """
    Convert a pyeda espresso boolean string into BNET format.
    Example: Or(And(~THY1, COL1A2), And(~THY1, ~CTSK))
    becomes: MYNODE, (!THY1 & COL1A2) | (!THY1 & !CTSK)
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

    # Build BNET line
    return f"{target}, {expr}"


def main():
    tree_text = """
|--- HLA-DMA <= 0.50
|   |--- HLA-DPA1 <= 0.50
|   |   |--- HLA-DQA1 <= 0.50
|   |   |   |--- CD74 <= 0.50
|   |   |   |   |--- class: 0
|   |   |   |--- CD74 >  0.50
|   |   |   |   |--- class: 0
|   |   |--- HLA-DQA1 >  0.50
|   |   |   |--- class: 0
|   |--- HLA-DPA1 >  0.50
|   |   |--- CD74 <= 0.50
|   |   |   |--- class: 1
|   |   |--- CD74 >  0.50
|   |   |   |--- class: 1
|--- HLA-DMA >  0.50
|   |--- CD74 <= 0.50
|   |   |--- HLA-DPA1 <= 0.50
|   |   |   |--- class: 0
|   |   |--- HLA-DPA1 >  0.50
|   |   |   |--- class: 1
|   |--- CD74 >  0.50
|   |   |--- HLA-DQA2 <= 0.50
|   |   |   |--- HLA-DPA1 <= 0.50
|   |   |   |   |--- class: 1
|   |   |   |--- HLA-DPA1 >  0.50
|   |   |   |   |--- HLA-DQA1 <= 0.50
|   |   |   |   |   |--- class: 1
|   |   |   |   |--- HLA-DQA1 >  0.50
|   |   |   |   |   |--- class: 1
|   |   |--- HLA-DQA2 >  0.50
|   |   |   |--- class: 1
    """.strip()

    # Optional: set the target gene name (could also parse from the filename or input)
    target_gene = "TARGET_GENE"

    # Step 1: Extract all class: 1 paths
    paths = extract_all_class1_paths(tree_text)

    # Step 2: Simplify Boolean logic
    simplified_rule = simplify_boolean_rule(paths)
    print(simplified_rule)

    # Step 3: Convert to BNET format
    bnet_rule = pyeda_to_bnet(simplified_rule, target=target_gene)
    print(bnet_rule)
    print("printed BNET rule")

if __name__ == "__main__":
    main()
