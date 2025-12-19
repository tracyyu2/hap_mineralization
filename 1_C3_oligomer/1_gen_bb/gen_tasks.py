import numpy as np
import os
from glob import glob
import random
import argparse
import pyrosetta
import json

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="A"/>
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
    <SecondaryStructureMetric name="ssmet" residue_selector="chainA"/>
    </SIMPLE_METRICS>
    <MOVERS>
    <MakePolyX name="makeA" aa="ALA" keep_pro="false"
        keep_gly="false" residue_selector="chainA" />
    </MOVERS>
    """)
    return obj


parser = argparse.ArgumentParser()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Path to input PDBs",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Path to write RFdiffusion outputs",
)

parser.add_argument(
    "--repeat_length_lower",
    default=25,
    type=int,
    help="Lower bound of repeat lengths to try",
)

parser.add_argument(
    "--repeat_length_upper",
    default=50,
    type=int,
    help="Upper bound of repeat lengths to try",
)

parser.add_argument(
    "--RFdiffusion_directory",
    default='/home/tracyy2/rf_diffusion_repeat/rf_diffusion_multichain',
    help="Path to RFdiffusion source",
)

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

fout = open('tasks','w')
cwd  = os.getcwd()

plist        = glob(args.input_directory+'/*pdb')
if len(plist)>250:
    plist_sample = random.sample(plist,250)
else:
    plist_sample = plist

pyrosetta.init("--mute all")
obj = initialize_xmls()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

for pdb in plist_sample:
    pose = pyrosetta.pose_from_file(pdb)
    ssmet = obj.get_simple_metric("ssmet")
    ss = ssmet.calculate(pose)
    
    rl = len(ss)/4 
    seedLoop = [1,rl]

    for rlT in range(args.repeat_length_lower,args.repeat_length_upper):
        ll = rlT-rl

        cstring="\\'A%d-%d,%d,A%d-%d,%d,A%d-%d,%d,A%d-%d,%d,0\\ B%d-%d,%d,B%d-%d,%d,B%d-%d,%d,B%d-%d,%d\\'"%(
                seedLoop[0],seedLoop[1],ll,seedLoop[0]+100,seedLoop[1]+100,ll,seedLoop[0]+200,seedLoop[1]+200,ll,seedLoop[0]+300,seedLoop[1]+300,ll,
                seedLoop[0],seedLoop[1],ll,seedLoop[0]+100,seedLoop[1]+100,ll,seedLoop[0]+200,seedLoop[1]+200,ll,seedLoop[0]+300,seedLoop[1]+300,ll)

        lline = "/software/containers/SE3nv.sif "+args.RFdiffusion_directory+"/run_inference.py "
        lline+= "contigmap.contigs=["+cstring+'] '
        lline+= "inference.output_prefix='"+args.output_directory+"/"+pdb.split('/')[-1][:-4]+"_%d"%rlT
        lline+= "' model.symmetrize_repeats=True model.repeat_length=%d model.symmsub_k=1 inference.num_designs=50 "%rlT
        lline+= "diffuser.T=50 inference.input_pdb='"+pdb+"' " 
        lline+= "inference.ckpt_override_path='"+args.RFdiffusion_directory+"/models/frozen_motif/BFF_5.pt' "
        lline+= "model.sym_method='max' model.copy_main_block_template=True model.main_block=1"

        fout.write(lline+'\n')

fout.close()
