import numpy as np
import pandas as pd
import pyrosetta as pr
import matplotlib.pyplot as plt

native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
occuring_experiments_fname = "/home/nieck/HSA_data/run_identifiers_everything"
input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
output_fname = "/home/nieck/HSA_data/results/matches_found"
fig_fname = "/home/nieck/HSA_data/results/true_matches_by_number_of_experiments.png"

#native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
#fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
#header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
#occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
#input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
#output_fname = "/home/niek/HSA_data/results/matches_found"
#fig_fname = "/home/niek/HSA_data/results/true_matches_by_number_of_experiments.png"


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

with open(output_fname, 'w') as out:
    out.write("number of experiments used and correct matches acquired by them\n")

matches = {}
for ex in runs:
    matches[ex] = set()

chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk.iterrows():
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
            matches[row['Run']].add((aa1Idx, aa2Idx))
    print("Finished chunk number %3d."%chunkCount)

idx = 0
allMatches = set()
for ex in matches:
    idx += 1
    allMatches.update(matches[ex])
    with open(output_fname, 'a') as out:
        out.write("%3d %4d\n"%(idx,len(allMatches)))

#mean number of matches per experiment: 869 (ex 1 2)

x = []
y = []
with open(output_fname, 'r') as result:
    for line in result:
        if line.startswith('number'):
            continue
        numbers = line.strip().split()
        x += [numbers[0]]
        y += [numbers[1]]

plt.figure()
plt.plot(x,y)
plt.xlabel("number of experiments used\n")
plt.ylabel("number of distinct correct matches found\n")
plt.savefig(fig_fname)
