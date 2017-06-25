import numpy as np
import pandas as pd
import pyrosetta as pr
import matplotlib.pyplot as plt
import contact_tools_fromMahmoud as ct
import crossval as xval
from sklearn.svm import SVC
import re

athome = False #use directory on my local computer
if not athome:
    native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/CompBioProject/headers/header_different_scoes'
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced_different_scores.csv"
    output_fname = "/home/nieck/HSA_data/results/crossval"
else:
    native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_different_scores'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_diffScores.csv"
    output_fname = "/home/niek/HSA_data/results/diffScores_classifier"

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
scorecolumns = ['betaCount', 'mgcScore', 'mgxScore', 'Precoursor Absolute Error', 'FragmentLibraryScore', 'AllScore', 'MatchScore', 'match score']
#['betaCount', 'mgcScore', 'mgxScore', 'Precoursor Absolute Error', 'MeanSquareError', 'FragmentLibraryScore', 'AllScore', 'MatchScore', 'match score']

#contacts = ct.extract_contacts_from_structure(native_fname)
#print("contacts:", contacts)

#take a number of rows out of each experiment
nofSamplesPerChunk = 70#10000
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
            X[pairkey] = np.zeros(3+len(scorecolumns))[np.newaxis,:]
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        X[pairkey][0,0] += 1 #number of contributing rows
        X[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
        for idx, scorekey in enumerate(scorecolumns):
            X[pairkey][0,idx+2] += row[scorekey]
        X[pairkey][0,-1] = res1pos.distance(res2pos) <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for pairkey in X:
    for idx in range(len(scorecolumns)):
        X[pairkey][0,idx+2] /= X[pairkey][0,0]
#for contact in contacts:
#    if (contact[0]+4, contact[1]+4) in X:
#        X[(contact[0]+4, contact[1]+4)][0,3] = True #contact

X = np.concatenate([X[x] for x in X], axis=0)
nofPositives = int(np.sum(X[:,-1]))
X = np.concatenate([X[X[:,-1]==1], X[X[:,-1]==0][:nofPositives]], axis=0)
print("training set shape:", X.shape)
print("proportion of positives in training set:", X[:,-1].mean())

#crossvalidate and train
classifier = xval.cv(X[:,:-1], X[:,-1], SVC, {'kernel':['linear', 'rbf']}, nfolds=5, nrepetitions=2, loss_function=xval.zero_one_loss)#xval.false_discovery_rate)
print('classifier kerneltype:', classifier.kernel)
print('xval loss (0-1):', classifier.cvloss)

#get a test set
Xtest = dict()
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
        if not pairkey in Xtest:
            Xtest[pairkey] = np.zeros(3+len(scorecolumns))[np.newaxis,:]
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        Xtest[pairkey][0,0] += 1 #number of contributing rows
        Xtest[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
        for idx, scorekey in enumerate(scorecolumns):
            Xtest[pairkey][0,idx+2] += row[scorekey]
        Xtest[pairkey][0,-1] = res1pos.distance(res2pos) <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for pairkey in Xtest:
    for idx in range(len(scorecolumns)):
        Xtest[pairkey][0,idx+2] /= Xtest[pairkey][0,0]
#for contact in contacts:
#    if (contact[0]+4, contact[1]+4) in Xtest:
#        Xtest[(contact[0]+4, contact[1]+4)][0,3] = True #contact
        
Xtest = np.concatenate([Xtest[x] for x in Xtest], axis=0)

print("test set shape:", Xtest.shape)
testPredictions = classifier.predict(Xtest[:,:-1])
print("proportion of positives in test set:", Xtest[:,-1].mean())
print("test set discoveries:", testPredictions.sum())
testNofPosPredicted = testPredictions.sum()
testNofFalseDiscoveries = testPredictions[np.logical_and(testPredictions==1, Xtest[:,-1]==0)].sum()
print("test set FDR:", float(testNofFalseDiscoveries)/testNofPosPredicted)
with open(output_fname, 'w') as f:
    f.write("classifier type: SVM, kernel: %s\n"%(classifier.kernel))
    f.write("test set discoveries: %d\n"%int(testPredictions.sum()))
    f.write("cvloss: %f\n"%classifier.cvloss)
