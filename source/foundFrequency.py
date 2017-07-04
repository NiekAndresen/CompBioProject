import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

athome = True
if athome:
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/HSA_data/header_reduced'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
    distance_fname = "/home/niek/HSA_data/1ao6/1ao6A.distances"
    output_fname = '/home/niek/HSA_data/results/foundFrequency'
else:
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/HSA_data/header_reduced'
    input_fname = "/home/nieck/HSA_data/data_experiment_1_2_reduced.csv"
    distance_fname = '/home/nieck/HSA_data/1ao6/1ao6A.distances'
    output_fname = '/home/nieck/HSA_data/results/foundFrequency'

# For a peptide link map gives three lists:
# - list of amino acids (1-letter) that might also be matches
# - list of probabilities that these aas are the match
# - list of index positions in the protein
def get_neighborhood_list(linkMap, startIdx):
    matches = re.findall(r'(\w+(?:-\w+)?)\((\d\.\d{2})\)', linkMap)
    aaList = []; valList = []; posList = []; origPosList = np.arange(len(matches)); valSum = 0
    for i,match in enumerate(matches):
        value = float(match[1])
        if len(match[0])==1:
            aaList += [match[0]]
            valList += [value]
            posList += [origPosList[i] + startIdx]
        elif match[0][-1].isupper():
            aaList += [match[0][-1]]
            valList += [value]
            posList += [origPosList[i] + startIdx]
        valSum += value
    valList = np.array(valList) / valSum
    return aaList, valList, posList

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

##checks several properties that a pair of indices needs to have in order to
##be taken into account
def valid_idx_pair(idx1, idx2):
    if idx1 < 5 or idx2 < 5: #these AAs are not in the .pdb
        return False
    if idx1-4 > 578 or idx2-4 > 578: #these AAs are not in the .pdb
        return False
    #if abs(idx1-idx2) < 12: #too low sequence separation
    #    return False
    if aa1Idx == aa2Idx:
        return False
    return True

pairdict = {} #storing the number of occurences for every pair of linkage sites
#determine frequency of occurence
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
chunkCount = 0
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk.iterrows():
        aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
        aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
        for i,aa1Idx in enumerate(pos1List):
            for j,aa2Idx in enumerate(pos2List):
                if val1List[i] * val2List[j] < .04: #very unlikely
                    continue
                site1 = int(aa1Idx)
                site2 = int(aa2Idx)
                if not valid_idx_pair(site1, site2):
                    continue
                if site1 > site2:
                    site2,site1 = site1,site2
                pairkey = (site1, site2)
                if not pairkey in pairdict:
                    pairdict[pairkey] = [0, False]
                pairdict[pairkey][0] += 1
                if row['MatchRank'] < 5: #this match has been top 5 in SOME scan
                    pairdict[pairkey][1] = True
    print("Finished chunk number %3d."%chunkCount)

toberemoved = []
for pair in pairdict: #remove all matches that have never been top 5
    if pairdict[pair][1] == False:
        toberemoved += [pair]
for pair in toberemoved:
    del pairdict[pair]

nofOccurrences = [pairv[0] for pairv in pairdict.values()]
labels = [dist[(pair[0]-4, pair[1]-4)] for pair in pairdict]
#print(nofOccurrences, labels)

with open(output_fname, 'w') as out:
    for i in range(len(labels)):
        out.write("%d, %f\n"%(nofOccurrences[i], labels[i]))
