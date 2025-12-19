import os
import glob

def combine_pdb_with_hetatm(protein_pdb_directory, mineral_pdb_file, output_directory):
    if not os.path.exists(mineral_pdb_file):
        print(f"Mineral PDB file not found: {mineral_pdb_file}")
        return

    with open(mineral_pdb_file, 'r') as min_file:
        mineral_lines = min_file.readlines()

    os.makedirs(output_directory, exist_ok=True)

    protein_pdb_files = glob.glob(os.path.join(protein_pdb_directory, '*.pdb'))

    if not protein_pdb_files:
        print(f"No PDB files found in the directory: {protein_pdb_directory}")
        return

    for protein_pdb_file in protein_pdb_files:
        if not os.path.exists(protein_pdb_file):
            print(f"Protein PDB file not found: {protein_pdb_file}")
            continue

        with open(protein_pdb_file, 'r') as prot_file:
            protein_lines = prot_file.readlines()

        file_name = os.path.basename(protein_pdb_file)
        combined_file_path = os.path.join(output_directory, f"combined_{file_name}")

        with open(combined_file_path, 'w') as combined_file:
            combined_file.writelines(protein_lines)
            combined_file.writelines(mineral_lines)
        
        print(f"Combined PDB file saved to: {combined_file_path}")

if __name__ == "__main__":
    protein_pdb_directory = "/home/tracyy2/git/tip_atom_mineralization/chr_diffusion_multichain/C3/4_ligandmpnn/010_d5_tube_PO4/input2"
    mineral_pdb_file = "/home/tracyy2/git/tip_atom_mineralization/chr_diffusion_multichain/C3/4_ligandmpnn/010_d5_tube_PO4/input2/010_d5_tube_PO4_shiftdown.pdb"
    output_directory = "/home/tracyy2/git/tip_atom_mineralization/chr_diffusion_multichain/C3/4_ligandmpnn/010_d5_tube_PO4/outputs"
    
    combine_pdb_with_hetatm(protein_pdb_directory, mineral_pdb_file, output_directory)
