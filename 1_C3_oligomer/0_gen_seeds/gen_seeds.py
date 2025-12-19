#backbone realignment and tranformation
#Alignment code adapted from: https://hunterheidenreich.com/posts/kabsch_algorithm/
import numpy as np
import itertools
from glob import glob
import shutil
import math
import os  # Import os for path operations
import pyrosetta
import argparse
import json

parser = argparse.ArgumentParser()

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="A"/>
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
    <SecondaryStructureMetric name="ssmet" residue_selector="chainA"/>
    <PerResidueClashMetric name="pclash" />
    </SIMPLE_METRICS>
    <MOVERS>
    <MakePolyX name="makeA" aa="ALA" keep_pro="false"
        keep_gly="false" residue_selector="chainA" />
    <DetectSymmetry name="detect" />
    </MOVERS>
    <FILTERS>
    <ClashCheck name="clasher" clash_dist="3.5"
        cutoff="4" 
         confidence="1.0" />
    </FILTERS>
    """)
    return obj

def remove_clashers(pdbname,obj,output_dir):
    for pdb_file in glob(output_dir+'/'+pdbname.split('/')[-1][:-4]+'*.pdb'):
        pose = pyrosetta.pose_from_file(pdb_file)
        pclash = obj.get_simple_metric("pclash")
        #clasher = obj.get_filter("clasher") ### MOD make a function that returns your residue selector
        #symm=obj.get_mover("detect")
        #symm.apply(pose)
        if sum(pclash.calculate(pose).values())>0:
            #print('removing '+pdb_file)
            os.remove(pdb_file)

def crystal_clash(pdbname,radius,output_dir):
    for pdb_in in glob(output_dir+'/'+pdbname.split('/')[-1][:-4]+'*.pdb'):
        fin = open(pdb_in,'r')
        in_line = fin.readline()
        while len(in_line)>0:
            if len(in_line)>4 and in_line[:4]=='ATOM':
                if in_line[13:16] in ['N  ','CA ','O  ','C  ']:
                    test_crds=np.zeros(2)
                    test_crds[0]=float(in_line[30:38])
                    test_crds[1]=float(in_line[38:46])
                    if np.sum(np.square(test_crds))<np.square(radius+4.3):
             #           print(pdb_in+' backbone clashes with crystal - removing')
                        os.remove(pdb_in)
                        fin.close()
                        break
            in_line= fin.readline()
        fin.close()

def align_points(A, B): 
    #This function aligns two arbirtrary sets of points
    # INPUTS:  A - first array of 3D points (N,3)
    #          B - second array of 3D points (N,3)
    # OUTPUTS: R - Rotation matrix to align  A to B (3,3)
    #          t - x,y,z translation to align A to B (1,3)
    #          rmsd - root mean square deviation of two points after alignment

    cen_A = np.mean(A, axis=0)
    cen_B = np.mean(B, axis=0)

    t = cen_A - cen_B

    a = A - cen_A
    b = B - cen_B

    H = np.dot(a.T, b)

    U, S, Vt = np.linalg.svd(H)

    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0

    R = np.dot(Vt.T, U.T)

    rmsd = np.sqrt(np.sum(np.square(np.dot(a, R.T) + cen_A - t - B)) / B.shape[0])

    return R, t, rmsd

def gen_motifs(R,t,xyz,cen_fit,twist,rise):
    # This function orients repeat motifs given alignment to crystal and superhelical parametwers
    # INPUTS:  R       - rotation matrix output by align points (3,3)
    #          t       - translation vector output by align points (1,3)
    #          xyz     - coordinate of atom to be aligned and superhelically transformed (1,3)
    #          cen_fit - centroid of full coordinate set used by align points (1,3)
    #          twist   - superhelical twist of repeats along the crystal
    #          rise    - superhelical rise of repeats along the crystal
    # OUTPUTS: xyzt    - transformed coordinates to output to motif pdb

    xyzt  = xyz - cen_fit #centers the coordinates
    xyzt  = np.dot(xyzt,R.T) + cen_fit - t #apply rotation and translation

    twist = twist*np.pi/180.0

    Rz    = np.array([[np.cos(twist),np.sin(twist),0],[-np.sin(twist),np.cos(twist),0],[0.0,0.0,1.0]])

    xyzt[:,2] = xyzt[:,2]+rise
    xyzt      = np.transpose(np.matmul(Rz,np.transpose(xyzt)))

    return xyzt

def retrieve_coords(pdb_in,ref_crystal):
    fit_crds = np.zeros((8,3))
    fin=open(pdb_in,'r')
    in_line = fin.readline()
    cnt=-2
    binder_inds=np.zeros(4)
    while len(in_line)>0:
        if len(in_line)>4 and in_line[:4]=='ATOM':
            if in_line[17:20] in ['GLU','ASP']:
                if ((in_line[17:20]=='GLU' and in_line[13:16]=='CD ') or
                    (in_line[17:20]=='ASP' and in_line[13:16]=='CG ')):
                    cnt+=2
                    fit_crds[cnt,0]=float(in_line[30:38])
                    fit_crds[cnt,1]=float(in_line[38:46])
                    fit_crds[cnt,2]=float(in_line[46:54])
                    binder_inds[int(cnt/2)]=int(in_line[22:26])
                elif (('d3' in pdb_in and in_line[13:16] in ['OE2','OD2']) or 
                      ('d5' in pdb_in and in_line[13:16] in ['OE1','OD1'])):
                    fit_crds[cnt+1,0]=float(in_line[30:38])
                    fit_crds[cnt+1,1]=float(in_line[38:46])
                    fit_crds[cnt+1,2]=float(in_line[46:54])
        in_line=fin.readline()
    fin.close()
 
    motif_bounds=[np.min(binder_inds)-3,np.max(binder_inds)+3]

    ref_crds = np.zeros((8,3))

    fin = open(ref_crystal,'r')
    in_line = fin.readline()
    i=-2
    while len(in_line)>0:
        if len(in_line)>4 and in_line[:4]=='ATOM':
            #query=int(in_line[22:26])
            #if query in ref_inds:
                if in_line[13:16] in ['CD ','CG ']:
                    i+=2
                    ref_crds[i,0]=float(in_line[30:38])
                    ref_crds[i,1]=float(in_line[38:46])
                    ref_crds[i,2]=float(in_line[46:54])
                elif (('d3' in pdb_in and in_line[13:16] in ['OE2','OD2']) or 
                      ('d5' in pdb_in and in_line[13:16] in ['OE1','OD1'])):
                    ref_crds[i+1,0]=float(in_line[30:38])
                    ref_crds[i+1,1]=float(in_line[38:46])
                    ref_crds[i+1,2]=float(in_line[46:54])
        in_line=fin.readline()
    fin.close()
    
    return fit_crds,ref_crds,motif_bounds

def gen_seeds(pdb_in,twists,rises,ref_crystal,outdir):
    # This function utilizes the above functions to align to crystal, then generate repeat chains
    # INPUTS:   pdb_in       - tip atom generated motif pdb
    #           twists,rises - lists of superhelical twists and rises to propagate repeats
    #           ref_cryal    - reference pdb of crystal used to align tip atom motifs
    #           outdir       - output directory to write the generated backbone motif pdbs

    (fit_crds, ref_crds, motif_bounds) = retrieve_coords(pdb_in,ref_crystal) 

    cen_fit = np.mean(fit_crds, axis=0)
    perms = itertools.permutations([0,1,2,3])

    for ea in perms:
        dum_crds = np.zeros((8,3)) 
        dum_crds[0,:] = ref_crds[ea[0]*2,:]
        dum_crds[2,:] = ref_crds[ea[1]*2,:]
        dum_crds[4,:] = ref_crds[ea[2]*2,:]
        dum_crds[6,:] = ref_crds[ea[3]*2,:]
        dum_crds[1,:] = ref_crds[ea[0]*2+1,:]
        dum_crds[3,:] = ref_crds[ea[1]*2+1,:]
        dum_crds[5,:] = ref_crds[ea[2]*2+1,:]
        dum_crds[7,:] = ref_crds[ea[3]*2+1,:]

        (R,t,rmsd)=align_points(fit_crds,dum_crds)

        if rmsd<0.5:
            break
            
    if rmsd>0.5:
        print('failed, motif might be bad: ',pdb_in)
        shutil.copyfile(pdb_in,'failures/'+pdb_in.split('/')[-1])
        return
    else:
        print(pdb_in, "succesfully aligned")

    for mm in range(len(rises)):
        rise=rises[mm]
        twist=twists[mm]

        pdb_out=pdb_in.split('/')[-1][:-4]+'_%4.2f_%4.2f.pdb'%(rise,twist)
        fout=open(outdir+'/'+pdb_out,'w')

        chains=['A','B','C','D','E','F','G','H']
        ctot=args.ctot

        for csym in range(args.ctot):
            for ll in range(0,4):
                resi=ll*100
                fin=open(pdb_in,'r')
                lline=fin.readline()
                xxx=np.zeros((1,3))
                while len(lline)>0:
                    if len(lline)>3 and lline[:4]=='ATOM' and int(lline[22:26])<100 and int(lline[22:26])>=(motif_bounds[0]) and int(lline[22:26])<=(motif_bounds[1]):
                        if lline[12:15]==' N ':
                            resi+=1
                        xxx[0,0]=float(lline[30:38])
                        xxx[0,1]=float(lline[38:46])
                        xxx[0,2]=float(lline[46:54])
                        twist_rep=(csym*(360.0/ctot)-twist*ll)
                        xxxN=gen_motifs(R,t,xxx,cen_fit,twist_rep,rise*ll)
                        fout.write(lline[:21]+chains[csym]+'%4d'%resi+lline[26:30]+'%8.3f%8.3f%8.3f'%(xxxN[0,0],xxxN[0,1],xxxN[0,2])+lline[54:])
                    lline=fin.readline()
                fin.close()

        fout.close()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Paths to PDBs for seed generation",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Paths to write output backbone seed PDBs",
)

parser.add_argument(
    "--reference_crystal",
    default='./',
    help="Path to PDB reference crystal used in aligment",
)

parser.add_argument(
    "--crystal_size",
    default='d3',
    help="crystal size (e.g. d3,d5,d7)",
)

parser.add_argument(
    "--ctot",
    type=int, #convert to integer
    default='6',
    help="protein symmetry",
)

args = parser.parse_args()


pyrosetta.init("-mute all")
obj = initialize_xmls()

# Define the output directory
output_dir = args.output_directory #'/net/scratch/tracyy2/git/0_seeds/07_09_2024_run1'
#output_dir = './'

# Create the directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)


#plist = glob('/net/scratch/tracyy2/HAp_binder_design_outputs/07_13_2024_design1/diffusion_out/010_d3_tube_OH_run1/*pdb')
plist = glob(args.input_directory+'/HAp*pdb')
#plist_sample=random.sample(plist,10)

for pdbname in plist:

    rise0  = 3.439
    twist0 = 60

    rises  = np.array([4*rise0,-4*rise0,rise0,3*rise0,rise0,3*rise0,-rise0,-3*rise0,-rise0,-3*rise0])
    twists = np.array([0,0,60,60,-60,-60.0,60,60,-60,-60.0])

    pdb_base = pdbname.split('/')[-1]
    ref_base = ''

    for i in pdb_base.split('_')[3:12]:
        ref_base=ref_base+i+'_'

    ref_base=ref_base[:-1]+'.pdb'

    ref_base_dir = ''

    for i in pdb_base.split('_')[3:7]:
        ref_base_dir=ref_base_dir+i+'_'

    ref_base_dir=ref_base_dir[:-1]+'/'

    ref_crystal=args.reference_crystal+'/'+ref_base_dir+ref_base

    if '_tube_OH' in pdbname:
        R0 = 8.2
    elif '_tube_Ca' in pdbname:
        R0 = 6.1
    elif '_tube_PO4' in pdbname:
        R0 = 10.2
    else:
        R0 = 0.0
    if '_d5_' in pdbname:
        R0 += 9.42

    gen_seeds(pdbname,twists,rises,ref_crystal,output_dir)
    remove_clashers(pdbname,obj,output_dir)
    crystal_clash(pdbname,R0,output_dir)
