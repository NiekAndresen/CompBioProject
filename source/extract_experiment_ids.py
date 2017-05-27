import pandas as pd

with open('/home/niek/HSA_data/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)
runs = set()
for chunk in chunks:
    for i,row in chunk.iterrows():
        runs.add(row['Run'])
nofruns = len(runs)
with open("/home/niek/HSA_data/Run_identifiers_1_2", 'w') as f:
    f.write(','.join(runs) + '\n')
