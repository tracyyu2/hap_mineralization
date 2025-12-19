import os

def rename_files(directory):
    # Check if the provided directory exists
    if not os.path.isdir(directory):
        print(f"The directory {directory} does not exist.")
        return
    
    # List all files in the directory
    for filename in os.listdir(directory):
        # Check if the filename starts with 'combined_'
        if filename.startswith("combined_"):
            # Create the new filename by removing the 'combined_' prefix
            new_filename = filename[len("combined_"):]

            # Construct full file paths
            old_file_path = os.path.join(directory, filename)
            new_file_path = os.path.join(directory, new_filename)
            
            # Rename the file
            os.rename(old_file_path, new_file_path)
            print(f"Renamed: {filename} -> {new_filename}")

# Specify the directory containing the files
directory = "/home/tracyy2/git/tip_atom_mineralization/chr_diffusion_multichain/C3/4_ligandmpnn/010_d5_tube_OH/outputs/unfixed_outputs/seqs/"

# Call the function to rename files
rename_files(directory)


