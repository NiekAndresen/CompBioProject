import numpy as np
import pandas as pd

athome = False #use directory on my local computer
if not athome:
    header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
else:
    header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"

with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

X = dict()
chunkCount = 0
count = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk.iterrows():
        if row['MatchRank'] > 5:
            continue
        else:
            count += 1
    print("Finished chunk number %3d."%chunkCount)
print("number of match rank <= 5: %d"%count)
