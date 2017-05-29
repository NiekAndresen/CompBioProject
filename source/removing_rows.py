
import pandas as pd
import numpy as np

with open('/home/niek/HSA_data/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#write header
with open("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", 'w') as f:
    f.write(','.join(columns) + '\n')

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2.csv", usecols=columns, chunksize=1e5)

for chunk in chunks:
    chunk[chunk['decoy']==0].to_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", mode='a', header=False, index=False)
    #ex_1_2 total rows: 800k
    #       decoy == 1: 600k; correlation: 0.0017
    #       decoy == 0: 200k; correlation: -0.2

#with open("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", 'r') as f:
#    header = f.readline()
#    print(header)
#    row = f.readline()
#    print(row)

