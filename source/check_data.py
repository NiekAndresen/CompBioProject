import numpy as np
import pandas as pd
import pyrosetta as pr

athome = False #use directory on my local computer
if not athome:
    native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
else:
    native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"

pr.init()
#load native protein
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()
native = pr.pose_from_pdb(native_fname)

with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

X = set()
chunkCount = 0
count = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk[chunk['MatchRank'] < 6].iterrows():
            aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
            aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
            for i, aa1Idx, aa2Idx in zip(range(len(pos1List)), pos1List, pos2List):
                aa1Idx = int(aa1Idx)
                aa2Idx = int(aa2Idx)
                if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                    continue
                if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                    continue
                if aa1Idx>aa2Idx:
                    aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
                pairkey = (aa1Idx, aa2Idx)
                if not pairkey in X:
                    X.add(pairkey)
    print("Finished chunk number %3d."%chunkCount)
print("number of different pairs occuring with rank <= 5: %d"%count)
