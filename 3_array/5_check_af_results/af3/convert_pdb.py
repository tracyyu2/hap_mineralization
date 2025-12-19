#!/usr/bin/env python3
import os
import glob
import shutil
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_file, pdb_file):
    """
    Convert a CIF file to a PDB file using BioPython.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def main():
    # Define directories
    outputs_dir = "/home/tracyy2/array_design/04_mlfold3/20251001"
    pdb_dir = os.path.join(outputs_dir, "pdb")
    
    # Create pdb directory if it doesn't exist
    os.makedirs(pdb_dir, exist_ok=True)

    # Find all CIF files in subdirectories of outputs
    cif_pattern = os.path.join(outputs_dir, "*", "*.cif")
    cif_files = glob.glob(cif_pattern)
    
    # Find all summary_confidences.json files in subdirectories of outputs
    json_pattern = os.path.join(outputs_dir, "*", "*summary_confidences.json")
    json_files = glob.glob(json_pattern)

    # Process CIF files: convert to pdb and save in pdb_dir
    for cif in cif_files:
        base_name = os.path.basename(cif)
        pdb_filename = os.path.splitext(base_name)[0] + ".pdb"
        pdb_path = os.path.join(pdb_dir, pdb_filename)
        try:
            convert_cif_to_pdb(cif, pdb_path)
            print(f"Converted {cif} to {pdb_path}")
        except Exception as e:
            print(f"Failed to convert {cif}: {e}")

    # Process JSON files: copy them to pdb_dir instead of moving
    for json in json_files:
        dest = os.path.join(pdb_dir, os.path.basename(json))
        try:
            shutil.copy2(json, dest)
            print(f"Copied {json} to {dest}")
        except Exception as e:
            print(f"Failed to copy {json}: {e}")

if __name__ == "__main__":
    main()
