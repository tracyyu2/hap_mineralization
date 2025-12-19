import numpy as np
from scipy.optimize import leastsq,minimize
from scipy.optimize import brute,basinhopping
from os import path
import os
from glob import glob
import math
import pyrosetta
import argparse
import json

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

def filter_clash(pdb_root,obj):
    pdb_file = glob(pdb_root+'*pdb')[0]
        
    pose    = pyrosetta.pose_from_file(pdb_file)
    clasher = obj.get_filter("clasher") ### MOD make a function that returns your residue selector

    symm    = obj.get_mover("detect")
    symm.apply(pose)

    if not clasher.apply(pose):
        print('removing '+pdb_file+' clashes')
        for pdb_file2 in glob(pdb_root+'*pdb'):
            os.remove(pdb_file2)
        return

    contacts = obj.get_filter("contacts")
    polyA = obj.get_mover("makeA")
    polyA.apply(pose)

    if contacts.score(pose)<10:
        print('removing '+pdb_file+' no contacts')
        for pdb_file2 in glob(pdb_root+'*pdb'):
            os.remove(pdb_file2)
        return

    breaks = obj.get_filter("cBreaks")

    if not breaks.apply(pose):
        print('removing '+pdb_file2+' chain breaks')
        for pdb_file2 in glob(pdb_root+'*pdb'):
            os.remove(pdb_file2)
        return


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

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="A"/>
            <Chain name="chainABC" chains="A,B,C"/>
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
    <SecondaryStructureMetric name="ssmet" residue_selector="chainA"/>
    </SIMPLE_METRICS>
    <MOVERS>
    <MakePolyX name="makeA" aa="ALA" keep_pro="false"
        keep_gly="false" residue_selector="chainABC" />
    <DetectSymmetry name="detect" />
    </MOVERS>
    <FILTERS>
    <ClashCheck name="clasher" clash_dist="3.5"
        cutoff="4"
         confidence="1.0" />
    <AtomicContactCount name="contacts" partition="chain"  distance="6.5"/>
    <ChainBreak name="cBreaks" chain_num="1" tolerance="0.3" threshold="0" />
    </FILTERS>
    """)
    return obj

def regen_repeats(tw,ri,xyz):
    Rz=np.array([[np.cos(-np.pi*tw/180.0),np.sin(-np.pi*tw/180),0],[-np.sin(-np.pi*tw/180),np.cos(-np.pi*tw/180),0],[0.0,0.0,1.0]])
    xyz[:,2]=xyz[:,2]+ri
    xyz=np.transpose(np.matmul(Rz,np.transpose(xyz)))
    return xyz

def align_motifs(R,t,xyz):
    xyzt=np.copy(xyz)
    cen_fit=np.mean(xyz, axis=0)
    xyzt=xyzt - cen_fit # #centers the coordinates
    xyzt=np.dot(xyzt,R.T) + cen_fit - t #apply rotation and translation
    return xyzt

def get_breaks(ss):
    Hstarts = []
    Hends   = []

    startR  = rl + 1
    startRt = startR

    for i in range(1,rl):
        if ss[startR-i] in 'HE':
            startRt-=1
        else:
            break

    startR = startRt

    in_helix = False 
    cnt      = 0 

    for i in range(startR,startR+rl):
        if ss[i] in 'HE':
            cnt+=1
        elif in_helix:
            Hends.append(i-1)
            in_helix = False
            cnt      = 0

        if cnt==3:
            Hstarts.append(i-2)
            in_helix = True

    if in_helix:
        Hends.append(startR+rl-1)

    if len(Hstarts)!=len(Hends):
        print("ss code failed to identify helices!")
        return

    init_size = Hends[-1]-Hstarts[0]+1
    SEs = []

    for i in Hstarts:
        for j in Hends:
            if (j-i+1)<init_size:
                SEs.append([i,j+rl])
            else:
                SEs.append([i,j])

    return SEs

def realign(pdbname,xyz,twist,rise,ref_crds,out_dir):
    splPDB = pdbname.split('/')

    (R,t,rmsd) = align_points(xyz,ref_crds)

    if rmsd>0.2:
        print(pdbname+' failed alignment',rmsd)
        return

    pose  = pyrosetta.pose_from_file(pdbname)
    ssmet = obj.get_simple_metric("ssmet")
    ss    = ssmet.calculate(pose)
    print(pdbname)
    SEs = get_breaks(ss)
    splPDB = pdbname.split('/')
    chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','O','P','Q','R','S']
    for bk in SEs:
      startR = bk[0]
      endR   = bk[1]
      mm     = 3
      tot    = float(mm)
      pdbOut = out_dir+'/'+splPDB[-1][:-4]+'_s%de%d_%dmer.pdb'%(startR,endR,mm)
      fout   = open(pdbOut,'w')
      resi   = 0

      for csym in range(mm):
          for ll in range(4):
              fin   = open(pdbname,'r')
              lline = fin.readline()
              xxx   = np.zeros((1,3))

              while len(lline)>0:
                   if len(lline)>3 and lline[:4]=='ATOM':
                       if ll<3 and (int(lline[22:26])<startR or int(lline[22:26])>=(startR+rl)):
                           lline = fin.readline()
                           continue

                       if ll==3 and (int(lline[22:26])<startR or int(lline[22:26])>(endR)):
                           lline = fin.readline()
                           continue

                       if lline[12:15]==' N ':
                            resi += 1

                       xxx[0,0] = float(lline[30:38])
                       xxx[0,1] = float(lline[38:46])
                       xxx[0,2] = float(lline[46:54])

                       xxxN0 = align_motifs(R,t,xxx)
                       xxxN  = regen_repeats(twist*ll+csym*360.0/tot,rise*ll,xxxN0,)

                       fout.write(lline[:21]+chains[csym]+'%4d'%resi+lline[26:30]+'%8.3f%8.3f%8.3f'%(xxxN[0,0],xxxN[0,1],xxxN[0,2])+lline[54:])

                   lline=fin.readline()

              fin.close()

      fout.close()

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Path to input PDBs",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Path to write output PDBs",
)

parser.add_argument(
    "--seed_directory",
    default='./',
    help="Path to reference seed PDBs",
)

parser.add_argument(
    "--task_id",
    default='./',
    help="Path to reference seed PDBs",
)

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

with open(args.output_directory+'/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

pyrosetta.init("--mute all")
obj = initialize_xmls()

l = glob(args.input_directory+'/*.pdb')

sublist =  l[int(args.task_id)::50]

for pdbname in sublist:
    xyz    = np.zeros((32,3))
    splPDB = pdbname.split('/')

    rl=int(pdbname.split('_')[-2])

    fin=open(pdbname,'r')
    i=0
    lline=fin.readline()

    while len(lline)>0:
        if len(lline)>3 and lline[:4]=='ATOM':
            if lline[17:20] in ['GLU','ASP'] and lline[13:15]=='CA':# and lline[21]=='A':
                xyz[i,0] += float(lline[30:38])
                xyz[i,1] += float(lline[38:46])
                xyz[i,2] += float(lline[46:54])

                i+=1

                if i>31:
                    break

        lline = fin.readline()

    fin.close()
    
    ref_pdb = ''

    for i in range(len(splPDB[-1].split('_'))-2):
        ref_pdb += splPDB[-1].split('_')[i]+'_'

    ref_pdb=args.seed_directory+'/'+ref_pdb[:-1]+'.pdb'

    ref_crds = np.zeros((32,3))
    fin      = open(ref_pdb,'r')

    i     = 0
    lline = fin.readline()

    while len(lline)>0:
        if len(lline)>3 and lline[:4]=='ATOM':
            if lline[17:20] in ['GLU','ASP'] and lline[13:15]=='CA':# and lline[21]=='A':
                ref_crds[i,0]+=float(lline[30:38])
                ref_crds[i,1]+=float(lline[38:46])
                ref_crds[i,2]+=float(lline[46:54])
                i+=1

                if i>31:
                    break

        lline=fin.readline()

    fin.close()
    rise=float(pdbname.split('_')[-4])
    twist=float(pdbname.split('_')[-3])

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

    realign(pdbname,xyz,twist,rise,ref_crds,args.output_directory)
    filter_clash(args.output_directory+'/'+pdbname.split('/')[-1][:-4],obj)
    crystal_clash(pdbname,R0,args.output_directory)
