
import pandas as pd
import numpy as np
import pyrosetta as pr
import matplotlib.pyplot as plt

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/niek/HSA_data/header_reduced'
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

pairdict = {} #storing the number of occurences for every pair of linkage sites

#determine frequency of occurence
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    for i,row in chunk.iterrows():
        site1 = int(row['ProteinLink1'])
        site2 = int(row['ProteinLink2'])
        pairkey = (site1, site2)
        if not pairkey in pairdict:
            pairdict[pairkey] = 0
        pairdict[pairkey] += 1

maximum = 0
for pair in pairdict:
    if pairdict[pair] > maximum:
        maximum = pairdict[pair]
print("maximum:")
print(maximum)

#determine chance to be true in crystal
frequencies = np.arange(maximum+1)
hitChance = np.zeros(maximum+1) #chance to be true in the crystal
nofPairsWFreq = np.zeros(maximum+1) #how many different pairs have this frequency of occurence
for pair in pairdict:
    nofPairsWFreq[pairdict[pair]] += 1
    if pair[0] < 5 or pair[1] < 5: #these AAs are not in the .pdb
        continue
    if pair[0]-4 > native.total_residue() or pair[1]-4 > native.total_residue():
        continue
    res1pos = native.residue(pair[0]-4).nbr_atom_xyz()
    res2pos = native.residue(pair[1]-4).nbr_atom_xyz()
    realdist = res1pos.distance(res2pos)
    hit = realdist > 15 and realdist < 25
    hitChance[pairdict[pair]] += hit
hitChance /= nofPairsWFreq

plt.figure()
plt.plot(frequencies, nofPairsWFreq)
plt.xlabel("appearance frequency")
plt.ylabel("number of different pairs")
plt.show()

antiHitChance = np.zeros(len(hitChance)) #to mark where there's no data
for i in range(len(hitChance)):
    if nofPairsWFreq[i] < 1:
        antiHitChance[i] = 1.
plt.figure()
plt.bar(frequencies, antiHitChance, color='gray', width=1, edgecolor="none")
plt.bar(frequencies, hitChance, color='blue', width=1, edgecolor="none")
plt.xlabel("appearance frequency")
plt.ylabel("chance to be true in the crystal")
plt.show()
