# coding: utf-8

import pandas as pd
import numpy as np
import pyrosetta as pr
import math

input_fname = "/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv"
native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A.fasta'
output_fname = "/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta"

pr.init()
#load native protein
with open(fasta_fname, 'r') as f:
    val = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            val += line.rstrip()
native_fasta = '_____' + val
native = pr.pose_from_pdb(native_fname)

with open(input_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
reconstruction = ['_'] * len(native_fasta) * 2
for chunk in chunks:
    for i,row in chunk.iterrows():
        try:
            if row['Linked AminoAcid 1'].isalpha() and len(row['Linked AminoAcid 1']) == 1:
                reconstruction[int(row['ProteinLink1'])] = row['Linked AminoAcid 1']
            if row['Linked AminoAcid 2'].isalpha() and len(row['Linked AminoAcid 2']) == 1:
                reconstruction[int(row['ProteinLink2'])] = row['Linked AminoAcid 2']
        except AttributeError:
            continue
        
reconstruction = ''.join(reconstruction)
reconstruction = reconstruction.strip('_')
with open(output_fname, 'w') as o:
    total_chars_written = 0
    while total_chars_written < len(reconstruction):
        if total_chars_written+80 >= len(reconstruction):
            o.write(reconstruction[total_chars_written:] + '\n')
        else:
            o.write(reconstruction[total_chars_written:total_chars_written+80] + '\n')
        total_chars_written += 80
