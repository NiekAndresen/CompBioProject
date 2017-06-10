import numpy as np
import pandas as pd
import pyrosetta as pr
import matplotlib.pyplot as plt

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/niek/HSA_data/header_reduced'
occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
output_fname = "/home/niek/HSA_data/results/using_truth"

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

with open(occuring_experiments_fname, 'r') as f:
    runs = f.readline()[:-1].split(',')
    
matches = {}
for ex in runs:
    matches[ex] = set()
    chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
    for chunk in chunks:
        for i,row in chunk.iterrows():
            if row['Run'] != ex:
                continue
            try:
                aa1Idx = int(row['ProteinLink1'])
                aa2Idx = int(row['ProteinLink2'])
            except ValueError:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                continue
            if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                continue
            res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
            res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
            realdist = res1pos.distance(res2pos)
            if row['MatchRank'] <= 5 and realdist < 25 and realdist > 15:
                matches[ex].add((aa1Idx, aa2Idx))

#mean number of matches per experiment: 869 (ex 1 2)

MatchesFoundWNofExperiments = np.zeros(len(matches))
newMatches = set(matches[list(matches.keys())[0]])
for i,key in enumerate(matches):
    if i==0: continue
    number = len(newMatches)
    newMatches.update(matches[key])
    MatchesFoundWNofExperiments[i] = MatchesFoundWNofExperiments[i-1] + len(newMatches) - number

plt.figure()
plt.plot(np.arange(len(matches)), MatchesFoundWNofExperiments)
plt.show()
