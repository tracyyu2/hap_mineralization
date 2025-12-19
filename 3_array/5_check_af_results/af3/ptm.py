import os
import json
import shutil
import re

def process_json_files(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop through all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".json"):
            input_file = os.path.join(input_dir, filename)
            
            # Open and read the JSON file
            with open(input_file, 'r') as f:
                data = json.load(f)
                
                # Check if the "ptm" value exists and if it's greater than 0.64
                ptm_value = data.get("ptm", 0)  # Default to 0 if ptm is not present
                
                # Extract the base name of the file (without path)
                base_name = filename.split('/')[-1]  # Get the filename from the full path
                
                # Remove "_summary_confidences" from the filename to generate name_stem
                name_stem = re.sub("_summary_confidences", '', base_name)
                name_stem = re.sub(".json", '', name_stem)
                name_stem = name_stem + "_model"
                
                # Print for debugging
                print(f"Base name: {base_name}")
                print(f"Name stem: {name_stem}")
                
                # Construct the corresponding pdb filename (correct extension)
                pdb_file = os.path.join(input_dir, f"{name_stem}.pdb")
                
                # Debugging output
                print(f"Looking for PDB file: {pdb_file}")
                
                # Check if ptm value is greater than setting, then copy the files
                if ptm_value > 0.67:
                    # Copy the JSON file to the output directory
                    shutil.copy(input_file, os.path.join(output_dir, filename))
                    shutil.copyfile(pdb_file, os.path.join(output_dir, f"{name_stem}.pdb"))
                    # Ensure the PDB file exists before attempting to copy it

# Example usage
input_dir = '/home/tracyy2/array_design/04_mlfold3/20251001/pdb'  # Replace with your input directory
output_dir = '/home/tracyy2/array_design/05_rmsd_iptm/20251001'  # Replace with your output directory

process_json_files(input_dir, output_dir)

