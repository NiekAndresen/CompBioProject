import numpy as np
import pandas as pd
import re

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

#load native protein
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

##checks several properties that a pair of indices needs to have in order to
##be taken into account
def valid_idx_pair(idx1, idx2):
    if idx1 < 5 or idx2 < 5: #these AAs are not in the .pdb
        return False
    if idx1-4 > 578 or idx2-4 > 578: #these AAs are not in the .pdb
        return False
    #if abs(idx1-idx2) < 12: #too low sequence separation
    #    return False
    if idx1 == idx2:
        return False
    return True

X = set()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk[chunk['MatchRank'] < 6].iterrows():
            aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
            aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
            for i, aa1Idx, aa2Idx in zip(range(len(pos1List)), pos1List, pos2List):
                if val1List[i] * val2List[i] < .1: #very unlikely
                continue
                aa1Idx = int(aa1Idx)
                aa2Idx = int(aa2Idx)
                if not valid_idx_pair(aa1Idx, aa2Idx):
                    continue
                if aa1Idx>aa2Idx:
                    aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
                pairkey = (aa1Idx, aa2Idx)
                if not pairkey in X:
                    X.add(pairkey)
    print("Finished chunk number %3d."%chunkCount)
print("number of unique pairs occuring with rank <= 5: %d"%len(X))
