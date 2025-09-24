import sys
import json

subprocess_path = sys.argv[0]
input_data_file_or_dir = sys.argv[1]
output_data_dir = sys.argv[2]
python_objects = sys.argv[3]


print(subprocess_path)
print(input_data_file_or_dir)
print(output_data_dir)

if len(python_objects) > 0:
    python_objects = python_objects.split("\n")
    for object in python_objects:
        python_objects[python_objects.index(object)] = json.loads(object.strip())

    print(python_objects)
    for object in python_objects:
        print(f"Object at index {python_objects.index(object)} of type {type(object)}: {object}")