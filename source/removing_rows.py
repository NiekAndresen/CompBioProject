
import pandas as pd
import numpy as np

athome = True
if not athome:
    header_input = '/home/nieck/CompBioProject/headers/header_reduced'
    input_file = "/home/nieck/HSA_data/SDA_HSA_Everything_3.csv"
    output_file = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
else:
    header_input = "/home/niek/Computational Biology/CompBioProject/headers/header_different_scores"
    input_file = "/home/niek/HSA_data/data_experiment_1_2.csv"
    output_file = "/home/niek/HSA_data/data_experiment_1_2_diffScores.csv"

#read in header to know what columns are supposed to be used
with open(header_input, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#write new csv:
#write header
with open(output_file, 'w') as f:
    f.write(','.join(columns) + '\n')
#write rows chunk wise
chunks = pd.read_csv(input_file, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunk = chunk[chunk['decoy']==0] #remove decoys
    chunk = chunk[lambda row: np.logical_and( [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 1']], [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 2']])] #remove self loops
    chunk.to_csv(output_file, mode='a', header=False, index=False)


