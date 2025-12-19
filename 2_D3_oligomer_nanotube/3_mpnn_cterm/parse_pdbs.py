from dateutil import parser
import numpy as np
import os, time, gzip, json
import glob 
import argparse
import pyrosetta
from pyrosetta.rosetta.core.select import get_residues_from_subset

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == "." and len(input_file) - i <= 5:
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

def get_residue_list(pose, selector):
    residue_list = get_residues_from_subset(selector.apply(pose))
    return residue_list

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Layer name="coreN2" select_core="true" select_boundary="false" select_surface="false" ball_radius="2.0" use_sidechain_neighbors="true" sc_neighbor_dist_exponent="0.7" core_cutoff="4.8"/>
            <Chain name="chainA" chains="A"/>
            <And name="cA" selectors="coreN2,chainA"/>
    </RESIDUE_SELECTORS>
    """)
    return obj

def make_pos_file(residue_list, pdb, out_dir, prefix = ""):
    pos_list = ''
    if '3mer' in pdb:
        divi=3
    elif '4mer' in pdb:
        divi=4
    elif '5mer' in pdb:
        divi=5
    elif '6mer' in pdb:
        divi=6

    llen=len(residue_list)
    spll=pdb.split('/')[-1].split('_')
    cnt=0
    for i in residue_list:
        pos_list += str(i) + "A,"
        cnt+=1
    filename = out_dir+'/pos/' + prefix + clean_file_name(pdb) + ".pos"
    pos_file = open(filename, "w")
    pos_file.write(pos_list)
    #print("Wrote " + filename)

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

argparser.add_argument("--input_directory", type=str, help="Path to a folder with pdb files, e.g. /home/my_pdbs/")
argparser.add_argument("--output_directory",   type=str, help="Path where to save .jsonl dictionary of parsed pdbs")

args = argparser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

if not os.path.exists(args.output_directory+'/pos'):
    os.makedirs(args.output_directory+'/pos')


with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

# folder_with_pdbs_path = args.pdb_folder
# save_path = args.out_path

#MODIFY THIS PATH TO YOUR FOLDER WITH PDBS
folder_with_pdbs_path = args.input_directory
#MODIFY OUTPUT PATH
save_path = args.output_directory + '/pdbs.jsonl'

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

aa_1_N = {a:n for n,a in enumerate(alpha_1)}
aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

def AA_to_N(x):
  # ["ARND"] -> [[0,1,2,3]]
  x = np.array(x);
  if x.ndim == 0: x = x[None]
  return [[aa_1_N.get(a, states-1) for a in y] for y in x]

def N_to_AA(x):
  # [[0,1,2,3]] -> ["ARND"]
  x = np.array(x);
  if x.ndim == 1: x = x[None]
  return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]


def parse_PDB_biounits(x, atoms=['N','CA','C'], chain=None):
  '''
  input:  x = PDB filename
          atoms = atoms to extract (optional)
  output: (length, atoms, coords=(x,y,z)), sequence
  '''
  xyz,seq,min_resn,max_resn = {},{},1e6,-1e6
  for line in open(x,"rb"):
    line = line.decode("utf-8","ignore").rstrip()

    if line[:6] == "HETATM" and line[17:17+3] == "MSE":
      line = line.replace("HETATM","ATOM  ")
      line = line.replace("MSE","MET")

    if line[:4] == "ATOM":
      ch = line[21:22]
      if ch == chain or chain is None:
        atom = line[12:12+4].strip()
        resi = line[17:17+3]
        resn = line[22:22+5].strip()
        x,y,z = [float(line[i:(i+8)]) for i in [30,38,46]]

        if resn[-1].isalpha(): 
            resa,resn = resn[-1],int(resn[:-1])-1
        else: 
            resa,resn = "",int(resn)-1
#         resn = int(resn)
        if resn < min_resn: 
            min_resn = resn
        if resn > max_resn: 
            max_resn = resn
        if resn not in xyz: 
            xyz[resn] = {}
        if resa not in xyz[resn]: 
            xyz[resn][resa] = {}
        if resn not in seq: 
            seq[resn] = {}
        if resa not in seq[resn]: 
            seq[resn][resa] = resi

        if atom not in xyz[resn][resa]:
          xyz[resn][resa][atom] = np.array([x,y,z])

  # convert to numpy arrays, fill in missing values
  seq_,xyz_ = [],[]
  try:
      for resn in range(min_resn,max_resn+1):
        if resn in seq:
          for k in sorted(seq[resn]): seq_.append(aa_3_N.get(seq[resn][k],20))
        else: seq_.append(20)
        if resn in xyz:
          for k in sorted(xyz[resn]):
            for atom in atoms:
              if atom in xyz[resn][k]: xyz_.append(xyz[resn][k][atom])
              else: xyz_.append(np.full(3,np.nan))
        else:
          for atom in atoms: xyz_.append(np.full(3,np.nan))
      return np.array(xyz_).reshape(-1,len(atoms),3), N_to_AA(np.array(seq_))
  except TypeError:
      return 'no_chain', 'no_chain'



pdb_dict_list = []

if folder_with_pdbs_path[-1]!='/':
    folder_with_pdbs_path = folder_with_pdbs_path+'/'


init_alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G','H', 'I', 'J','K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T','U', 'V','W','X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g','h', 'i', 'j','k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't','u', 'v','w','x', 'y', 'z']
extra_alphabet = [str(item) for item in list(np.arange(300))]
chain_alphabet = init_alphabet + extra_alphabet

biounit_names = glob.glob(folder_with_pdbs_path+'*.pdb')

def parser(biounit_names,save_path,init_alphabet,extra_alphabet,chain_alphabet):
    c=0

    cnt=0
    fin=open(save_path,'w')
    for biounit in biounit_names:
        cnt+=1
        if cnt%10==0:
            print(cnt)
        my_dict = {}
        s = 0
        concat_seq = ''
        concat_N = []
        concat_CA = []
        concat_C = []
        concat_O = []
        concat_mask = []
        coords_dict = {}
        for letter in chain_alphabet:
            xyz, seq = parse_PDB_biounits(biounit, atoms=['N','CA','C','O'], chain=letter)
            if type(xyz) != str:
                concat_seq += seq[0]
                my_dict['seq_chain_'+letter]=seq[0]
                coords_dict_chain = {}
                coords_dict_chain['N_chain_'+letter]=xyz[:,0,:].tolist()
                coords_dict_chain['CA_chain_'+letter]=xyz[:,1,:].tolist()
                coords_dict_chain['C_chain_'+letter]=xyz[:,2,:].tolist()
                coords_dict_chain['O_chain_'+letter]=xyz[:,3,:].tolist()
                my_dict['coords_chain_'+letter]=coords_dict_chain
                s += 1
        fi = biounit.rfind("/")
        my_dict['name']=biounit[(fi+1):-4]
        my_dict['num_of_chains'] = s
        my_dict['seq'] = concat_seq
        if s < len(chain_alphabet):
       # pdb_dict_list.append(my_dict)
            fin.write(json.dumps(my_dict) + '\n')
            c+=1
        
    fin.close()       
    print(f'Successfully finished: {len(biounit_names)} pdbs')

def chain(save_path):
    with open(save_path, 'r') as json_file:
        json_list = list(json_file)

    my_dict = {}
    for json_str in json_list:
        result = json.loads(json_str)
        all_chain_list = [item[-1:] for item in list(result) if item[:9]=='seq_chain'] #['A','B', 'C',...]
        masked_chain_list = all_chain_list
        visible_chain_list = []
        my_dict[result['name']]= (masked_chain_list, visible_chain_list)

    with open(args.output_directory+'/pdbs_masked.jsonl', 'w') as f:
        f.write(json.dumps(my_dict) + '\n')


def fix(save_path):
    with open(save_path, 'r') as json_file:
        json_list = list(json_file)

    my_dict = {}

    for json_str in json_list:
        result = json.loads(json_str)

        fixed_position_dict = {}

        # Process each PDB file in the input directory
        pdb_file_path = os.path.join(args.input_directory, f"{result['name']}.pdb")
        with open(pdb_file_path, 'r') as fin:
            chain_residues = {}  # To collect residues for each chain
            chain_a_residues = []  # To store first 61 residues of Chain A
            for lline in fin:
                if len(lline) > 3 and lline[:4] == 'ATOM' and lline[13:15] == 'CA':  # Ensure it's a CA atom
                    chain_id = lline[21]  # Get chain identifier
                    res_num = int(lline[22:26].strip())

                    # If we are reading Chain A, collect the first 61 residues
                    if chain_id == 'A' and len(chain_a_residues) < 61:
                        chain_a_residues.append(res_num)

            # After processing the PDB, apply the first 61 residues of Chain A to all chains
            # Overwrite residue numbers (apply the same first 61 residues to all chains)
            for chain_id in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':  # Assuming chains are labeled A, B, C, etc.
                fixed_position_dict[chain_id] = chain_a_residues

        my_dict[result['name']] = fixed_position_dict

    # Write the updated fixed positions to the output file
    with open(os.path.join(args.output_directory, 'pdbs_fixed.jsonl'), 'w') as f:
        f.write(json.dumps(my_dict) + '\n')


def get_sel():
    pyrosetta.init("--mute all")

    obj = initialize_xmls()

    for pdb_file in glob.glob(args.input_directory+'/*pdb'):
        pose = pyrosetta.pose_from_file(pdb_file)

        residue_selector = obj.get_residue_selector("cA") ### MOD make a function that returns your residue selector
        residue_list = get_residue_list(pose, residue_selector)

        make_pos_file(residue_list, pdb_file, args.output_directory)


import json

def make_tie(save_path):
    with open(save_path, 'r') as json_file:
        json_list = list(json_file)

    my_dict = {}
    for json_str in json_list:
        result = json.loads(json_str)
        all_chain_list = sorted([item[-1:] for item in result if item.startswith('seq_chain')])
        tied_positions_list = []

        # Read the positions from the .pos file
        with open(f"{args.output_directory}/pos/{result['name']}.pos") as fin:
            lline = fin.readline()
        
        # Split and parse positions
        spl = lline.split(',')[:-1]
        Tspl = spl

        # Create tied positions list based on parsed positions
        for i in spl:
            if i in Tspl:
                dupes = [int(j[:-1]) for j in Tspl if j[:-1] == i[:-1]]
                T2spl = [j for j in Tspl if j[:-1] != i[:-1]]
                temp_dict = {chain: dupes for chain in all_chain_list}
                tied_positions_list.append(temp_dict)
                Tspl = T2spl

        # Fill remaining positions for each chain
        chain_length = len(result["seq_chain_A"])
        splNum = [int(i[:-1]) for i in spl]
        for i in range(1, chain_length + 1):
            if i not in splNum:
                temp_dict = {chain: [i] for chain in all_chain_list}
                tied_positions_list.append(temp_dict)

        my_dict[result['name']] = tied_positions_list

    # Write final dictionary to output JSON file
    with open(args.output_directory+'/pdbs_tied.jsonl', 'w') as f:
        f.write(json.dumps(my_dict) + '\n')
 
parser(biounit_names,save_path,init_alphabet,extra_alphabet,chain_alphabet)
chain(save_path)
fix(save_path)
get_sel()
make_tie(save_path)


print('Finished')
