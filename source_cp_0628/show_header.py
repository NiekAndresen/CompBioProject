
import pandas as pd
import numpy as np

with open('/home/nieck/CompBioProject/headers/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#write header
#with open("/home/nieck/HSA_data/data_experiment_1_2_nodecoys.csv", 'w') as f:
#    f.write(','.join(columns) + '\n')

#chunks = pd.read_csv("/home/nieck/HSA_data/SDA_HSA_Everything_3.csv", usecols=columns, chunksize=1e5)

#for chunk in chunks:
#    chunk[chunk['decoy']==1].to_csv("/home/nieck/HSA_data/data_experiment_1_2_nodecoys.csv", mode='a', header=False, index=False)

with open("/home/nieck/HSA_data/sda_filtered_kolja.csv", 'r') as f:
    header = f.readline()
    print(header)
    row = f.readline()
    print(row)

