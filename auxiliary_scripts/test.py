import sys
import json
import helper_functions as hf

def import_cmd_args(constant_count):
    """
    Returns input data path, output data path and each element of a list of python objects received
    from `execute_subprocess` as command line arguments. Set `constant_count` to the number of
    variables you expect to receive as command line arguments.

    Example
    -------

    >>> input_data_path, output_data_path, annotations_list, verbose,  = import_cmd_args(4)
    """

    if len(constant_count) < 2:
        raise ValueError("Please provide at least 2 constant names to import_cmd_args. The first one is the input data path, the second one is the output data path.")

    # assign python objects list
    if len(constant_count) > 2 and sys.argv[3] != "":
        python_objects = sys.argv[3]
    elif len(constant_count) > 2 and sys.argv[3] == "":
        raise ValueError("Trying to assign python objects, but no python objects list was provided. Please provide a list with python objects to execute_subprocess.")
    elif len(constant_count) == 2 and sys.argv[3] != "":
        raise ValueError("Trying to assign python objects, but no constant_names for them where provided in the subprocess. Please provide constant_names for the expected python objects to impor_cmd_args.")

    # import python objects
    python_objects = python_objects.split("\n")
    for object in python_objects:
        python_objects[python_objects.index(object)] = json.loads(object.strip())

    return sys.argv[1], sys.argv[2], *python_objects # * unpacks the list into separate elements

