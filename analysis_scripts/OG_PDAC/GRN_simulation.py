# script should take as argument what modules to run
# simulate GRN (asynchronous or synchronous updates)
# attractor analysis
# KO analysis

# NOTE: large networks lead to stack overflow in NuSMV
# NOTE: shortening the network at this point is not possible, because then nodes that are used as inputs might not exist anymore
# NOTE: apparently it is a problem in pyboolnet that "." is not allowed in node names, so we encode it as _DOT_
# NOTE: tarjan's algorithm is only viable for very small networks, as requires and stg, so we must use random walks or minimal trap space approximations
# NOTE: random seeds for random walks from pyboolnet can also set constant 0 nodes to 1, so the first and second state per random walk can be impossible 
#       states. As long as and attractor is not reached within the first or second state, this does not matter though. I tried this by limiting walk length to 2 and no attractors where found, so it should be fine.

from pyboolnet import file_exchange, attractors, trap_spaces
import helper_functions as hf

# import cmd args 
input_data_file, output_data_dir, verbose = hf.import_cmd_args(3)
vprint = hf.make_vprint(verbose)

# open bnet file and sanitize "." in node names
vprint("Importing bnet...")
with open(input_data_file, "r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace(".", "_DOT_")   
bnet_string = "".join(lines)

# find prime implicants
vprint("Converting bnet to prime implicants...")
primes = file_exchange.bnet2primes(bnet_string)

def find_unique_attractors_rw(primes, update_method: str ="synchronous", walk_length:int =1000, walk_attempts: int =10, n_walks: int = 1, seed: dict = None) -> list[dict]:
    """
    Find unqiue attractors by randomly walking through state space. The random walk stops once
    an attractor is found. As such it only finds one attractor per walk and is not guaranteed
    to find all attractors in n walks. A higher n_walks will increase the chance of 
    finding more unqiue attractors. Will raise an error if no attractor is found in the given
    number of walk_attempts.

    Parameters:
        primes (list): list of prime implicants of network to search
        update_method (str): "synchronous" or "asynchronous"
        walk_length (int): max length of each random walk
        walk_attempts (int): max number of random walks to find each attractor 
        n_walks (int): number of walks to finish / non - unique attractors to find
        seed (dict): a state of the network to start the random walk from (only 1 walk is done if a seed is provided)

    Returns:
        unique_attractors (list): list of unique attractors, each attractor is a dictionary
        (key = node name, value = 0 or 1)
    """
    attractors_list = []

    if seed is None:
        # find attractors with random walks
        vprint(f"Finding {n_walks} attractor by random walk...")
        for i in range(n_walks):
            attractor_rw = attractors.find_attractor_state_by_randomwalk_and_ctl(primes, update=update_method, length=walk_length, attempts=walk_attempts) # raises an error if no attractor found in attempts tries
            attractors_list.append(attractor_rw)
    else:
        vprint("Finding attractors from seed...")
        attractor_rw = attractors.find_attractor_state_by_randomwalk_and_ctl(primes, update=update_method, length=walk_length, attempts=walk_attempts, initial_state=seed) # raises an error if no attractor found in attempts tries
        attractors_list.append(attractor_rw)

    # find unique attractors in list
    vprint("Identifying unique attractors...")
    unique_attractors = []
    for attractor in attractors_list:
        if attractor not in unique_attractors:
            unique_attractors.append(attractor)

    print(f"Found {len(unique_attractors)} unique attractors.")

    return unique_attractors

def find_unique_attractors_mts(primes, update_method: str ="synchronous", n_trapspaces: any = 10000) -> list[dict]:
    """
    Find unqiue attractors using minimal trap space approximation. Not all attractors
    are guaranteed to be in an MTS (so called motif avoidant attractors), an MTS may 
    contain multiple attractors, and an attractor in an MTS can be much smaller than
    the MTS. Thus 3 criteria exist for the quality of an MTS as an approximation for
    an attractor:     

    completeness: whether every attractor is contained in one of the network's minimal
    trap spaces

    univocality: whether each minimal trap spaces contains exactly one attractor

    faithfulness: whether the number of free variables of the minimal trap space is 
    equal to the number of oscillating variables of the attractors contained in it

    The function checks all minimal trap spaces for univocality and faithfulness. 
    If they are both true for a trap space, its attractor is added to a list.
    Unique attractors are then extracted from the list.

    Completeness is reported to stdout.

    Parameters:
        primes (list): list of prime implicants of network to search
        update_method (str): "synchronous" or "asynchronous"
        n_trapspaces (int or convertible to int): max number of minimal trap spaces to check.
            cannot be above 2.147e9 or will cause integer overflow.

    Returns:
        unique_attractors (list): list of unique attractors, each attractor is a dictionary
        (key = node name, value = 0 or 1)
    """
    n_trapspaces = int(n_trapspaces) # standardize format so input like 1e6 also works

    # find minimal trap spaces
    vprint("Finding minimal trap spaces...")
    min_t_spaces = trap_spaces.compute_trap_spaces(primes, max_output=n_trapspaces)
    attractors_list = []

    # check univocality and faithfulness for each minimal trap space
    vprint("Finding attractors in minimal trap spaces...")
    for trap_space in min_t_spaces:
        univocal = attractors.univocality(primes, update=update_method, trap_space=trap_space)
        faithful = attractors. faithfulness(primes, update=update_method, trap_space=trap_space)
        if univocal and faithful:
            attractors_list.append(trap_space)


    # check if minimal trap spaces contain all attractors
    vprint("Checking completeness...")
    complete = attractors.completeness(primes, update=update_method, max_output=n_trapspaces)

    if complete:
        print("All attractors found in minimal trap spaces.")
    else:
        print("Not all attractors found in minimal trap spaces.")
    
    # find unique attractors in list
    vprint("Identifying unique attractors...")
    unique_attractors = []
    for attractor in attractors_list:
        if attractor not in unique_attractors:
            unique_attractors.append(attractor)

    print(f"Found {len(unique_attractors)} unique attractors.")

    return unique_attractors

# print output (none, so executor does not cry about missing output)
print("Output: " + "none")


if __name__ == "__main__":

    seed = {'AC013275-DOT-1': 1, 'AC013275_DOT_1': 1, 'AC091182_DOT_1': 0, 'ADGRG1': 1, 'ADH1C': 1, 'AIF1L': 0, 'ALDOB': 1, 'AMY2B': 1, 'ANXA9': 1, 'APCS': 0, 'APOH': 1, 'AQP1': 1, 'ARC': 0, 'BICC1': 1, 'CA12': 0, 'CA2': 1, 'CAPN5': 1, 'CCR10': 0, 'CDH11': 0, 'CDH13': 0, 'CFB': 1, 'CLDN2': 1, 'CLEC14A': 0, 'CLTB': 1, 'CPPED1': 0, 'CRIP2': 1, 'CTSF': 1, 'CTSS': 1, 'CYSTM1': 1, 'DCTPP1': 1, 'DUSP26': 0, 'EHF': 1, 'ESRP1': 1, 'FLRT2': 1, 'FSTL3': 0, 'GALNT12': 1, 'GATA2': 0, 'GLIPR2': 0, 'HIST1H2AI': 0, 'HMCN1': 0, 'HSPG2': 0, 'HYDIN': 0, 'IFI27': 1, 'IFITM3': 1, 'IGFN1': 0, 'KISS1': 0, 'KRT19': 1, 'LARGE2': 0, 'LCP1': 0, 'LDB3': 0, 'LGALS7B': 0, 'LINC01480': 0, 'LOX': 0, 'LTBP2': 0, 'LTC4S': 1, 'MAD2L1': 0, 'MATK': 0, 'MEG8': 0, 'MGP': 0, 'MS4A4A': 0, 'MS4A8': 1, 'MUSTN1': 0, 'MXRA8': 0, 'MYH11': 0, 'MYOCD': 0, 'NCAM1': 0, 'NECTIN1': 0, 'NID1': 0, 'NRP2': 0, 'PCYOX1': 1, 'PDGFRA': 0, 'PGF': 0, 'PHGR1': 0, 'PLAAT3': 1, 'PLAT': 0, 'PMAIP1': 1, 'PODN': 0, 'PROCR': 1, 'PRSS21': 0, 'RRP12': 0, 'SCGB2A1': 1, 'SFRP2': 0, 'SHANK3': 0, 'SLC28A3': 1, 'SLC39A5': 1, 'SLIT3': 0, 'SLITRK5': 0, 'SOCS3': 1, 'SPRR1B': 0, 'ST18': 0, 'STARD10': 1, 'STRA6': 0, 'SYNPO2': 0, 'TACSTD2': 1, 'TBX2': 0, 'THSD7A': 0, 'TM4SF18': 1, 'TNFRSF1B': 0, 'WNT2': 0, 'WWC1': 1, 'Z99572_DOT_1': 0}

    try:
        attractors_mts = find_unique_attractors_mts(primes, update_method="synchronous")
        # attractors_mts = find_unique_attractors_rw(primes, update_method="synchronous", n_walks=100)
        for attractor in attractors_mts:
            print(attractor)
    except:
        raise RuntimeError("No attractors found")
 