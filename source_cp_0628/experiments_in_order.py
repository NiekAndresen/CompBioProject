import pandas as pd

input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"

with open(occuring_experiments_fname, 'r') as f:
    runs = f.readline()[:-1].split(',')
nof_experiments = len(runs)
print("number of experiments:", nof_experiments) #2

chunks = pd.read_csv(input_fname, chunksize=1e5)
nof_ex_changes = 0
currentEx = None
for chunk in chunks:
    for i,row in chunk.iterrows():
        if not currentEx:
            currentEx = row['Run']
        if currentEx != row['Run']:
            print("change from %s to %s."%(currentEx, row['Run']))
            nof_ex_changes += 1
            currentEx = row['Run']
print("number of experiment changes:", nof_ex_changes) #191

#RESULT: data is NOT sorted by experiment
