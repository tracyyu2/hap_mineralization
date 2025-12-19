import os
from glob import glob
import Bio.PDB
import pyrosetta
import argparse
import shutil
import json

simple_threading_mover = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover()
pyrosetta.init("--beta_nov16 true --mute all")

def threader(fs, templ_dir, out_dir):
    fin = open(fs, 'r')
    lline = fin.readline()
    lline = fin.readline()
    tes1 = lline.split(':')[0] * 6
    parent_pose = pyrosetta.pose_from_file(templ_dir + "/" + fs.split('/')[-1][:-5] + ".pdb")
    simple_threading_mover.set_sequence(tes1, 1)
    simple_threading_mover.apply(parent_pose)
    parent_pose.dump_pdb(out_dir + '/' + fs.split('/')[-1][:-3] + ".pdb")

parser = argparse.ArgumentParser()

parser.add_argument("--input_directory", default='./', help="Path to MPNN outputs")
parser.add_argument("--output_directory", default='./', help="Path to write colabfold outputs")
parser.add_argument("--reference_directory", default='./', help="Path to reference PDBs used for RMSD")

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

if not os.path.exists(args.output_directory + '/fastas'):
    os.makedirs(args.output_directory + '/fastas')

with open(args.output_directory + '/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

print(len(glob(args.input_directory + '/outputs*/*pdb')))

for pdb in sorted(glob(args.input_directory + '/outputs*/*pdb')):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    # Load structures
    ref_structure = pdb_parser.get_structure("reference", args.reference_directory + '/' + pdb[:-65].split('/')[-1] + '.pdb')
    sample_structure = pdb_parser.get_structure("sample", pdb)

    def get_ca_atoms(structure, chain_ids):
        ca_atoms = []
        for chain_id in chain_ids:
            chain = structure[0][chain_id]
            for res in chain:
                if 'CA' in res:
                    ca_atoms.append(res['CA'])
        return ca_atoms

    # Select chains A, B, C, etc, change
    ref_ca_atoms = get_ca_atoms(ref_structure, ['A'])
    sample_ca_atoms = get_ca_atoms(sample_structure, ['A'])

    # Calculate RMSD
    rmsd = 1000
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_ca_atoms, sample_ca_atoms)
    super_imposer.apply(sample_structure[0].get_atoms())
    rmsd = super_imposer.rms

    print(f"RMSD: {rmsd:.4f} for {pdb}")

    if rmsd < 3.0:
        in_ind = pdb.split('/')[-2][7:]
        fasta_file = args.input_directory + '/fastas' + in_ind + '/' + pdb[:-63].split('/')[-1] + '.fa'
        shutil.copyfile(fasta_file, args.output_directory + '/fastas/' + pdb[:-63].split('/')[-1] + '.fa')
        threader(fasta_file, args.reference_directory, args.output_directory)
        print(f"RMSD {rmsd:.4f} passed: {pdb}")

