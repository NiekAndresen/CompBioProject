import pandas as pd
import numpy as np

athome = False
if not athome:
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
    occuring_experiments_fname = "/home/nieck/HSA_data/run_identifiers_everything"
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
    output_fname = "/home/nieck/HSA_data/results/everything_statistics"
    distance_fname = '/home/nieck/HSA_data/1ao6/1ao6A.distances'
else:
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/HSA_data/header_reduced'
    occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
    output_fname = "/home/niek/HSA_data/results/ex_1_2_statistics"
    distance_fname = "/home/niek/HSA_data/1ao6/1ao6A.distances"

#load native fasta
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()

#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

threshold = 20.#Angstrom
rankstrue = np.zeros(11) #for each rank the number of times a row with this rank was true
ranknum = np.zeros(11) #for each rank the number of times a row with this rank occured
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
chunkCount = 0
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
        if aa1Idx-4 > 578 or aa2Idx-4 > 578:
            continue
        if aa1Idx==aa2Idx:
            continue
        if row['MatchRank'] > 10:
            rank = 11 - 1 #-1 for idx
        else:
            rank = row['MatchRank']
        ranknum[rank] += 1
        if aa1Idx > aa2Idx:
            aa1Idx,aa2Idx = aa2Idx,aa1Idx
        if dist[(aa1Idx-4, aa2Idx-4)] < 20.:
            rankstrue[rank] += 1
    print("Finished chunk number %3d."%chunkCount)

print(rankstrue / ranknum)


