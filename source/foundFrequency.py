import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

athome = True
if athome:
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/HSA_data/header_reduced'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
    distance_fname = "/home/niek/HSA_data/1ao6/1ao6A.distances'
else:
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/HSA_data/header_reduced'
    input_fname = "/home/nieck/HSA_data/data_experiment_1_2_reduced.csv"
    distance_fname = '/home/nieck/HSA_data/1ao6/1ao6A.distances'

#load native fasta
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()

with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

pairdict = {} #storing the number of occurences for every pair of linkage sites

#determine frequency of occurence
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    for i,row in chunk.iterrows():
        site1 = int(row['ProteinLink1'])
        site2 = int(row['ProteinLink2'])
        pairkey = (site1, site2)
        if not pairkey in pairdict:
            pairdict[pairkey] = [0, False]
        pairdict[pairkey][0] += 1
        if row['MatchRank'] < 5: #this match has been top 5 in SOME scan
            pairdict[pairkey][1] = True
        #    pairdict[pairkey] = 0
        #pairdict[pairkey] += 1

toberemoved = []
for pair in pairdict: #remove all matches that have never been top 5
    if pairdict[pair][1] == False:
        toberemoved += [pair]
for pair in toberemoved:
    del pairdict[pair]

maximum = 0
for pair in pairdict:
    if pairdict[pair][0] > maximum:
        maximum = pairdict[pair][0]
print("maximum:")
print(maximum)

#determine chance to be true in crystal
frequencies = np.arange(maximum+1)
hitChance = np.zeros(maximum+1) #chance to be true in the crystal
nofPairsWFreq = np.zeros(maximum+1) #how many different pairs have this frequency of occurence
for pair in pairdict:
    nofPairsWFreq[pairdict[pair][0]] += 1
    if pair[0] < 5 or pair[1] < 5: #these AAs are not in the .pdb
        continue
    if pair[0]-4 > 578 or pair[1]-4 > 578:
        continue
    hit = dist[(aa1Idx-4, aa2Idx-4)] <= 20
    hitChance[pairdict[pair][0]] += hit
hitChance /= nofPairsWFreq

#plt.figure()
#plt.bar(frequencies[50:], nofPairsWFreq[50:])
#plt.xlabel("appearance frequency")
#plt.ylabel("number of different pairs")
#plt.show()

antiHitChance = np.zeros(len(hitChance)) #to mark where there's no data
for i in range(len(hitChance)):
    if nofPairsWFreq[i] < 1:
        antiHitChance[i] = .2
plt.figure()
plt.bar(frequencies, antiHitChance, color='gray', width=1, edgecolor="none")
plt.bar(frequencies, hitChance, color='blue', width=1, edgecolor="none")
plt.xlabel("appearance frequency")
plt.ylabel("chance to be true in the crystal")
plt.show()
