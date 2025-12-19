
import os

def rename_pdb_files(directory):
    """
    Renames all PDB files in the given directory by adding a trailing underscore
    before the file extension.

    Args:
        directory (str): Path to the directory containing PDB files.
    """
    # Get a list of all PDB files in the directory
    pdb_files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    
    # Loop through each PDB file and rename it
    for pdb_file in pdb_files:
        # Construct the new filename with a trailing underscore
        new_name = pdb_file[:-4] + '_' + pdb_file[-4:]  # Removes '.pdb' and adds '_.' before '.pdb'
        old_path = os.path.join(directory, pdb_file)
        new_path = os.path.join(directory, new_name)
        
        # Rename the file
        os.rename(old_path, new_path)
        print(f"Renamed: {old_path} -> {new_path}")

# Example usage
if __name__ == "__main__":
    directory = '/home/tracyy2/git/tip_atom_mineralization/chr_diffusion_multichain/C3/3_realign/08_22_2024_d5_PO4'  
    # Replace with the path to your directory
    rename_pdb_files(directory)


