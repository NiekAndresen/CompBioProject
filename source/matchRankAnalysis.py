
import pandas as pd
import numpy as np
import pyrosetta as pr

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/niek/HSA_data/header_reduced'
occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
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

with open(occuring_experiments_fname, 'r') as f:
    runs = f.readline()[:-1].split(',')

for ex in runs:
    scans = {}
    matchRanks = {}
    chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
    for chunk in chunks:
        for i,row in chunk.iterrows():
            if row['Run'] != ex:
                continue
            if not row['Scan'] in scans.keys():
                scans[row['Scan']] = []
            try:
                aa1Idx = int(row['ProteinLink1'])
                aa2Idx = int(row['ProteinLink2'])
            except ValueError:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                continue
            if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                continue
            scans[row['Scan']] += [(row['match score'], row['MatchRank'])]#, row['Autovalidation'])]
    for scan in scans.keys():
        if not all(scans[scan][i][0] >= scans[scan][i+1][0] for i in range(len(scans[scan])-1)):#if not sorted
            scans[scan].sort(key=lambda match: match[0], reverse=True) #descending by score
    if len(scans) == 0:
        continue
    print("experiment:", ex)
    sortedbyscore = [False] * len(scans)
    sortedbyrank = [False] * len(scans)
    #autovalidated = [False] * len(scans)
    for i,scan in enumerate(scans.keys()):
        sortedbyscore[i] = all(scans[scan][j][0] >= scans[scan][j+1][0] for j in range(len(scans[scan])-1))
        sortedbyrank[i]  = all(scans[scan][j][1] <= scans[scan][j+1][1] for j in range(len(scans[scan])-1))
        #autovalidated[i] = any(scans[scan][j][2] for j in range(len(scans[scan])))
    print("all sorted by score:", all(sortedbyscore)) #True
    print("all sorted by rank:", all(sortedbyrank)) #True
    #print("all any autovalidated:", all(autovalidated))
    print("ratio scans sorted by rank:", sum(sortedbyrank)/len(sortedbyrank)) #1.0
    #print("ratio scans sorted by rank out of not any autovalidated:", sum([x and not y for (x,y) in zip(sortedbyrank, autovalidated)])/(len(scans)-sum(autovalidated))) #~50%

#RESULT: rank 1 means highest score, rank2 means second highest score etc. ranks can appear multiple times
#column 'Autovalidation' and 'validated' are still a mystery




