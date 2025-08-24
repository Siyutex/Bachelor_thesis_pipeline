import os
import requests
import time

# Path to the directory containing your downloaded folders
base_dir = r"C:\Users\Julian\Downloads"

# List all folders (assumes each folder name is a file UUID)
folder_uuids = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]

# Loop through each folder UUID
for folder_uuid in folder_uuids:
    api_url = f"https://api.gdc.cancer.gov/files/{folder_uuid}"
    
    try:
        # get the case_id for each folder UUID
        r = requests.get(api_url, params={"fields": "cases.case_id"})
        r.raise_for_status()
        case_id = r.json()['data']['cases'][0]['case_id']  # get the first case_id from the response
        print("case_id: " + case_id) # print the value of the case_id field from the json object the request returned 
        
        # get the disease type and primary site for each case_id
        r = requests.get(f"https://api.gdc.cancer.gov/cases/" + case_id, params={"fields": "disease_type,primary_site"})
        disease_type = r.json()['data']['disease_type']
        primary_site = r.json()['data']['primary_site']
        print("disease_type: " + disease_type)  # print the value of the disease_type field from the json object the request returned
        print("primary_site: " + primary_site)  # print the value of the primary

        # change folder name to include disease type and primary site
        new_folder_name = f"{disease_type}_{primary_site}_{folder_uuid}"
        os.rename(os.path.join(base_dir, folder_uuid), os.path.join(base_dir, new_folder_name))

        # Be polite to the API
        time.sleep(0.1)
        
    except Exception as e:
        print(f"Failed for {folder_uuid}: {e}")