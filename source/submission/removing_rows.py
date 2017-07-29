## Niek Andresen for Computational Biology Project Summer Term 2017

# This script reads through a big dataset from CLMS like
# 'SDA_HSA_Everything_3.csv' and removes self-loops, decoys as well as all
# columns that are not requested.
# It requires a file that only contains a header with column names that are
# supposed to be kept (e.g. 'header_reduced.csv').

# GIVE PATHS:
header_input = "/home/nieck/CompBioProject/headers/header_reduced.csv"
input_file = "/home/nieck/HSA_data/SDA_HSA_Everything_3.csv"
output_file = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"

import pandas as pd
import numpy as np

#read in header to know what columns are supposed to be used
with open(header_input, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#print(columns)

#write new csv:
#write header
with open(output_file, 'w') as f:
    f.write(','.join(columns) + '\n')
#write rows chunk wise
chunks = pd.read_csv(input_file, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunk = chunk[chunk['decoy']==0] #remove decoys
    list_valid_AA_1 = [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 1']] #self loops have longer descriptions, not just one AA letter
    list_valid_AA_2 = [str(e).isalpha() and len(str(e))==1 for e in row['Linked AminoAcid 2']]
    chunk = chunk[lambda row: np.logical_and(list_valid_AA_1, list_valid_AA_2)] #remove self loops
    chunk.to_csv(output_file, mode='a', header=False, index=False)
