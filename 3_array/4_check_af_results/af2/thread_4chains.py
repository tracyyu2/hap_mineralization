import os 
from glob import glob
import Bio.PDB
import pyrosetta
import argparse
import shutil
import json
simple_threading_mover = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover()
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta import *
from pyrosetta.rosetta import *
detect_symmetry = pyrosetta.rosetta.protocols.symmetry.DetectSymmetry()

def threader(fs,templ_dir,out_dir):
    fin=open(fs,'r')
    lline=fin.readline()
    lline=fin.readline()
    tes1=lline.split(':')[0]+lline.split(':')[0]+lline.split(':')[0]+lline.split(':')[0] #change if chain numbers change
    parent_pose = pyrosetta.pose_from_file(templ_dir+"/"+fs.split('/')[-1][:-5]+".pdb")
    #detect_symmetry.apply(parent_pose)
    simple_threading_mover.set_sequence(tes1,1)
    simple_threading_mover.apply(parent_pose)
    parent_pose.dump_pdb(out_dir+'/'+fs.split('/')[-1][:-3]+".pdb")
                                                        
pyrosetta.init("--beta_nov16 true --mute all")

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Path to MPNN outputs",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Path to write colabfold outputs",
)

parser.add_argument(
    "--reference_directory",
    default='./',
    help="Path to reference PDBs used for RMSD",
)

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

if not os.path.exists(args.output_directory+'/fastas'):
    os.makedirs(args.output_directory+'/fastas')

with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)


print(len(glob(args.input_directory+'/output*/*pdb')))
for pdb in sorted(glob(args.input_directory+'/outputs*/*pdb')):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)

    ref_structure = pdb_parser.get_structure("reference", pdb)
    sample_structure = pdb_parser.get_structure("sample", args.reference_directory+'/'+pdb[:-65].split('/')[-1]+'.pdb')

    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]

    ref_list = []

    for ref_chain in ref_model:
        ref_atoms_temp=[]
        for ref_res in ref_chain:
            ref_atoms_temp.append(ref_res['CA'])
        ref_list.append(ref_atoms_temp)
        
    sample_list = []

    for sample_chain in sample_model:
        sample_atoms_temp=[]
        for sample_res in sample_chain:
            sample_atoms_temp.append(sample_res['CA'])
        sample_list.append(sample_atoms_temp)

    rmsd=1000
    super_imposer = Bio.PDB.Superimposer()
    from itertools import combinations
    in_list = list(combinations(range(len(ref_list)), 2))
    for i in in_list:
        ref_atoms=ref_list[i[0]]+ref_list[i[1]]
        sample_atoms=sample_list[0]+sample_list[1]
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms()) 
        if rmsd>super_imposer.rms:
            rmsd=super_imposer.rms

    print(rmsd)
    if rmsd<6.0:
        in_ind = pdb.split('/')[-2][7:]
        shutil.copyfile(args.input_directory+'/fastas'+in_ind+'/'+pdb[:-63].split('/')[-1]+'.fa',args.output_directory+'/fastas/'+pdb[:-63].split('/')[-1]+'.fa')
        threader(args.input_directory+'/fastas'+in_ind+'/'+pdb[:-63].split('/')[-1]+'.fa',args.reference_directory,args.output_directory)
        print('%6.4f'%rmsd,pdb)

