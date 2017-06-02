
import pandas as pd
import numpy as np
import pyrosetta as pr

pr.init()
#load native protein
native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
with open('/home/niek/HSA_data/1ao6/1ao6A.fasta', 'r') as f:
    val = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            val += line.rstrip()
native_fasta = '_____' + val
native = pr.pose_from_pdb(native_fname)

with open('/home/niek/HSA_data/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)
for chunk in chunks:
    for i,row in chunk.iterrows():
        try:
            aa1Idx = int(row['ProteinLink1'])
            aa2Idx = int(row['ProteinLink2'])
            aa1Name = row['Linked AminoAcid 1']
            aa2Name = row['Linked AminoAcid 2']
            if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                continue
            if aa1Idx-4 < 1 or aa2Idx-4 < 1:
                continue
            if aa1Idx-28 > native.total_residue() or aa2Idx-28 > native.total_residue():
                continue
            if aa1Idx-28 < 1 or aa2Idx-28 < 1:
                continue
            print('offset of 5:')
            print(aa1Name, native.residue_type(aa1Idx-4).name1(), native.sequence()[aa1Idx-5], native_fasta[aa1Idx])
            #print('offset of 28:')
            #print(aa1Name, native.residue_type(aa1Idx-27).name1(), native.sequence()[aa1Idx-28], native_fasta[aa1Idx-23])
            #print(native.residue_type(aa2Idx-4).name1(), native.sequence()[aa2Idx-5], aa2Name, native_fasta[aa2Idx])
            #print("in fasta offset 28:", native_fasta[int(row['ProteinLink1'])-28])
            #print("in fasta offset 5:", native_fasta[int(row['ProteinLink1'])-5])
            #print("in csv:", row['Linked AminoAcid 1'])
        except ValueError:
            continue

#with open('/home/niek/HSA_data/1ao6/1ao6A_seq_from_pdb.fasta', 'w') as f:
#    f.write(native.sequence() + '\n')

#probability of having at least one pair of residues in contact if two peptides appear as row
#count = 0
#success = 0
#chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)
#for chunk in chunks:
#    for i,row in chunk.iterrows():
#        if(i>1000): break
#        foundOne = False
#        indexPairWarValid = False
#        count += 1
#        try:
#            for aa1Idx in range(int(row['Start1'])-4,int(row['Start1'])-4+int(row['LengthPeptide1'])):
#                if foundOne: break
#                for aa2Idx in range(int(row['Start2'])-4,int(row['Start2'])-4+int(row['LengthPeptide2'])):
#                    if aa1Idx < 1 or aa2Idx < 1: #what AAs are those? I don't have them in the .pdb and .fasta
#                        continue
#                    if aa1Idx > native.total_residue() or aa2Idx > native.total_residue():
#                        continue
#                    indexPairWasValid = True
#                    if foundOne: break
#                    res1pos = native.residue(aa1Idx).nbr_atom_xyz()
#                    res2pos = native.residue(aa2Idx).nbr_atom_xyz()
#                    realdist = res1pos.distance(res2pos)
#                    if realdist < 20.:
#                        success += 1
#                        foundOne = True
#        except ValueError:
#            count -= 1
#            continue
#        if not indexPairWasValid:
#            count -= 1
#print(float(success)/count)



