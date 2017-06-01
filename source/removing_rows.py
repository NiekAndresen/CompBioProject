
import pandas as pd
import numpy as np

header_input = '/home/niek/HSA_data/header_reduced'
input_file = "/home/niek/HSA_data/data_experiment_1_2.csv"
output_file = "/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv"

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
    chunk[chunk['decoy']==0].to_csv(output_file, mode='a', header=False, index=False)
    #todo remove self loops


