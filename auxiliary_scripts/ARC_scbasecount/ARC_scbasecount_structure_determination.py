import gcsfs
import json

basepath = "gs://arc-scbasecount/2025-02-25/"
fs = gcsfs.GCSFileSystem(project=basepath)

def map_file_system(path):
    result = {}

    for entry in fs.ls(path=path):
        full_path = entry
        print ("entry: " + entry)
        print("path: " + path)

        # "." to check if it's a file or folder, "/" count to limit the depth of the recursion, and avoid infinite recursion by checking if the entry is the same as the current path
        if "." not in full_path and path.count("/") < 6 and entry != "arc-scbasecount/2025-02-25/" and (entry + "/") != path:
            print(full_path + "/")
            result[entry] = map_file_system(full_path + "/")
        if "." in full_path:
            result[entry] = None
    return result

with open("folder_structure_2.json", "w") as f:
    json.dump(map_file_system(basepath), f, indent=3)

# folder_structure.json has "/" count < 5 and folder_structure_2.json has "/" count < 6
# folder_structure_2.json is a bit longer so there are some extra folders that weren't caught