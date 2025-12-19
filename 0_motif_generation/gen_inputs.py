from itertools import permutations

def write_motif_pdb(dN,surf,tips,ag):

    res_dict = {'D':'ASP' , 'E':'GLU'}
    c_dict   = {'D':'CG  ', 'E':'CD  '}
    o1_dict  = {'D':'OD1 ', 'E':'OE1 '}
    o2_dict  = {'D':'OD2 ', 'E':'OE2 '}

    fout     = open('input_motif/010_'+dN+'_tube_'+surf+
                     '_%d_%d_%d_%d'%(tips[0],tips[1],tips[2],tips[3])+'_'+ag+'.pdb','w')
    resid    = 0

    fin     = open('/home/tracyy2/from/Harley/tubes/010_'+dN+'_tube_'+surf+'_tips2.pdb','r')
    in_line = fin.readline()

    while len(in_line)>0:
        if len(in_line)>4 and in_line[:4]=='ATOM' and int(in_line[22:26]) in tips:
            if in_line[13:15] in ['CD','CG']:
                resid+=1
                dummy_bb = ("ATOM      1  N   %s B%4d      17.430  39.052  10.587  1.00 14.60\n"%(res_dict[ag[resid-1]],resid)+
                            "ATOM      2  CA  %s B%4d      17.326  40.045  11.654  1.00 12.60\n"%(res_dict[ag[resid-1]],resid)+
                            "ATOM      3  C   %s B%4d      18.590  40.114  12.523  1.00 14.86\n"%(res_dict[ag[resid-1]],resid)+
                            "ATOM      4  O   %s B%4d      18.641  40.865  13.504  1.00 15.76\n"%(res_dict[ag[resid-1]],resid))
                fout.write(dummy_bb)
                fout.write(in_line[:13]+c_dict[ag[resid-1]]+res_dict[ag[resid-1]]+in_line[20:22]+'%4d'%resid+in_line[26:])
            elif in_line[13:16] in ['OE1','OD1']:
                fout.write(in_line[:13]+o1_dict[ag[resid-1]]+res_dict[ag[resid-1]]+in_line[20:22]+'%4d'%resid+in_line[26:])
            elif in_line[13:16] in ['OE2','OD2']:
                fout.write(in_line[:13]+o2_dict[ag[resid-1]]+res_dict[ag[resid-1]]+in_line[20:22]+'%4d'%resid+in_line[26:])
        in_line=fin.readline()
    fin.close()

    fin     = open('/home/tracyy2//from/Harley/tubes/010_'+dN+'_tube_'+surf+'.pdb','r')
    in_line = fin.readline()

    while len(in_line)>0:
        fout.write(in_line)
        in_line=fin.readline()
    fin.close()

    fout.close()


DEperms = ['DDDD','EDDD','DEDD','EEDD','DDED','EDED','DEED','EEED','DDDE','EDDE','DEDE','EEDE','DDEE','EDEE','DEEE','EEEE']
for i in DEperms:
    write_motif_pdb('d5','PO4',[28,563,164,162],i)
for i in DEperms:
    write_motif_pdb('d5','PO4',[28,563,478,704],i)
#for i in DEperms:
    #write_motif_pdb('d5','OH',[141,205,48,47],i)
#for i in DEperms:
    #write_motif_pdb('d5','Ca',[276,273,347,343],i)





