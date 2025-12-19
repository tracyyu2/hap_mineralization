import pyrosetta
import shutil
import glob
import argparse
import os 
import json

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="A"/>
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
    <SecondaryStructureMetric name="ssmet" residue_selector="chainA"/>
    </SIMPLE_METRICS>
    <FILTERS>
      <ChainBreak name="cBreaks" chain_num="1" tolerance="0.3" threshold="5" />
    </FILTERS>
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

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

pyrosetta.init("-mute all")
obj = initialize_xmls()

cnt = 0

for pdb_file in glob.glob(args.input_directory+'/*pdb'):
    print(pdb_file)
    cnt+=1

    if cnt%100==0:
        print(cnt)
    pose = pyrosetta.pose_from_file(pdb_file)
 
    cBreak = obj.get_filter("cBreaks")
    ssmet  = obj.get_simple_metric("ssmet")

    if cBreak.apply(pose):
        ss = ssmet.calculate(pose)

        if float(ss.count('H')+ss.count('E'))/len(ss)>0.75:
            shutil.copyfile(pdb_file,args.output_directory+'/'+pdb_file.split('/')[-1])

