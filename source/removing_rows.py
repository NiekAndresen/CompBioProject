
import pandas as pd
import numpy as np

header_input = '/home/nieck/CompBioProject/headers/header_reduced'
input_file = "/home/nieck/HSA_data/SDA_HSA_Everything_3.csv"
output_file = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"

with open(header_input, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#write header
with open(output_file, 'w') as f:
    f.write(','.join(columns) + '\n')

chunks = pd.read_csv(input_file, usecols=columns, chunksize=1e5)

for chunk in chunks:
    chunk = chunk[chunk['decoy']==0]
    chunk = chunk[lambda row: np.logical_and( [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 1']], [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 2']])]
    chunk.to_csv(output_file, mode='a', header=False, index=False)


