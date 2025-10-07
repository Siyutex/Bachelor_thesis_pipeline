# script should take as argument what modules to run
# simulate GRN (asynchronous or synchronous updates)
# attractor analysis
# KO analysis

# NOTE: shortening the network at this point is not possible, because then nodes that are used as inputs might not exist anymore
# NOTE: apparently it is a problem in pyboolnet that "." is not allowed in node names, so we encode
from pyboolnet import file_exchange, attractors
import helper_functions as hf

# import cmd args 
input_data_file, output_data_dir, verbose = hf.import_cmd_args(3)
vprint = hf.make_vprint(verbose)


# sanitize "." in node names
with open(input_data_file, "r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace(".", "_DOT_")
    
bnet_string = "".join(lines)
# write to txt
with open("gode.txt", "w") as f:
    f.write(bnet_string)

# import GRN
primes = file_exchange.bnet2primes(bnet_string)

attractor = attractors.find_attractor_state_by_randomwalk_and_ctl(primes, update="asynchronous", length=1000)
print(attractor)