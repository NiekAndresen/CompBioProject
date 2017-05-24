
# coding: utf-8

# In[56]:

import pandas as pd
import numpy as np
import pyrosetta as pr


# In[58]:

pr.init()
#load native protein
native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
with open('/home/niek/HSA_data/1ao6/1ao6A.fasta', 'r') as f:
    key = ''
    for line in f:
        if line.startswith('>'):
            if key:
                break
            key, val = line[1:].rstrip().split()[0], ''
        elif key:
            val += line.rstrip()
native_fasta = val
native = pr.pose_from_pdb(native_fname)


# In[24]:

with open('/home/niek/HSA_data/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)


# In[55]:

#i = 0
#for chunk in chunks:
#    if i > 0:
#        break
#    i += 1
#    print(chunk.iloc[0])


# In[50]:

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2.csv", usecols=columns, chunksize=1e5)


# In[47]:

for chunk in chunks:
    for i,row in chunk.iterrows():
        if i > 0:
            break
        if row['decoy'] == 1:
            continue
        aa1Idx = int(row['ProteinLink1'])
        aa2Idx = int(row['ProteinLink2'])
        aa1Name = row['Linked AminoAcid 1']
        aa2Name = row['Linked AminoAcid 2']
        print(row['BasePeptide1'], native_fasta[int(row['Start1']-5):int(row['Start1']-5)+int(row['LengthPeptide1'])])
        print(row['BasePeptide2'], native_fasta[int(row['Start2']-5):int(row['Start2']-5)+int(row['LengthPeptide2'])])
        print(native_fasta[aa1Idx-5], aa1Name, row['BasePeptide1'])
        print(native_fasta[aa2Idx-5], aa2Name, row['BasePeptide2'])


# In[32]:




# In[ ]:



