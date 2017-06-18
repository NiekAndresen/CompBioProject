import numpy as np
import pandas as pd
import pyrosetta as pr
import matplotlib.pyplot as plt
import contact_tools_fromMahmoud as ct
import crossval as xval
from sklearn.svm import SVC

native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
occuring_experiments_fname = "/home/nieck/HSA_data/run_identifiers_everything"
input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
output_fname = "/home/nieck/HSA_data/results/crossval"

#native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
#fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
#header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
#occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
#input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
#output_fname = "/home/niek/HSA_data/results/crossval"


pr.init()
#load native protein
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()
native = pr.pose_from_pdb(native_fname)

with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

with open(occuring_experiments_fname, 'r') as f:
    runs = f.readline()[:-1].split(',')

contacts = ct.extract_contacts_from_structure(native_fname)
#print("contacts:", contacts)

#take a number of rows out of each experiment
nofSamplesPerChunk = 25#10000
X = dict()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    rows = chunk.sample(n=nofSamplesPerChunk)
    for i,row in rows.iterrows():
        try:
            aa1Idx = int(row['ProteinLink1'])
            aa2Idx = int(row['ProteinLink2'])
        except ValueError:
            continue
        if aa1Idx > aa2Idx:
            aa1Idx,aa2Idx = aa2Idx,aa1Idx#smaller index first for contact list
        if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
            continue
        if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
            continue
        pairkey = (aa1Idx, aa2Idx)
        if not pairkey in X:
            X[pairkey] = np.zeros(5)[np.newaxis,:]
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        X[pairkey][0,0] += 1#number of contributing rows
        X[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
        X[pairkey][0,2] += row['match score'] #score
        X[pairkey][0,-1] = res1pos.distance(res2pos) <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for contact in contacts:
    if (contact[0]+4, contact[1]+4) in X:
        X[(contact[0]+4, contact[1]+4)][0,3] = True #contact

X = np.concatenate([X[x] for x in X], axis=0)
print(X.shape)
print(X[:,-1].mean())

classifier = xval.cv(X[:,:-1], X[:,-1], SVC, {'kernel':['linear']}, nfolds=5, nrepetitions=2)
print('xval loss:', classifier.cvloss)
with open(output_fname, 'w') as f:
    f.write("%f\n"%classifier.cvloss)
