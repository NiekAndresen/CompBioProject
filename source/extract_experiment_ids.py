import pandas as pd

header_fname = '/home/niek/HSA_data/header_reduced'
input_fname = "/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv"
output_fname = "/home/niek/HSA_data/run_identifiers_1_2"

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
