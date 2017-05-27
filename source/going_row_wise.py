# coding: utf-8

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

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2.csv", usecols=columns, chunksize=1e5)

for chunk in chunks:
    for i,row in chunk.iterrows():
        if i > 0:
            break
        if row['decoy'] == 0:
            continue
        aa1Idx = int(row['ProteinLink1']) #index in the fasta string from above (preceding '_____')
        aa2Idx = int(row['ProteinLink2'])
        aa1Name = row['Linked AminoAcid 1']
        aa2Name = row['Linked AminoAcid 2']
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        print(res1pos.distance(res2pos))
        
        
        
        
        #print(native.residue_type(aa1Idx-4).name1(), native.sequence()[aa1Idx-5], aa1Name, native_fasta[aa1Idx])
        #print(native.residue_type(aa2Idx-4).name1(), native.sequence()[aa2Idx-5], aa2Name, native_fasta[aa2Idx])
        #above: because pose numbering of residues starts at 1 instead of 0
        
        #print(native.total_residue(), len(native.sequence()), len(native_fasta))
        #print(row['BasePeptide1'], native_fasta[int(row['Start1']):int(row['Start1'])+int(row['LengthPeptide1'])])
        #print(row['BasePeptide2'], native_fasta[int(row['Start2']):int(row['Start2'])+int(row['LengthPeptide2'])])
        #print(native_fasta[aa1Idx], aa1Name, row['BasePeptide1'])
        #print(native_fasta[aa2Idx], aa2Name, row['BasePeptide2'])



