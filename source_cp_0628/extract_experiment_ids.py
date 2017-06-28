import pandas as pd

header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_3.csv"
output_fname = "/home/nieck/HSA_data/run_identifiers_everything"

with open(header_fname, 'r') as hr:
    header = hr.readline().rstrip() #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
runs = set()
for chunk in chunks:
    for i,row in chunk.iterrows():
        runs.add(row['Run'])
nofruns = len(runs)
with open(output_fname, 'w') as f:
    f.write(','.join(runs) + '\n')
