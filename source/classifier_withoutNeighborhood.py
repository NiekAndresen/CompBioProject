import numpy as np
import pandas as pd
import pyrosetta as pr
import matplotlib.pyplot as plt
import contact_tools_fromMahmoud as ct
import crossval as xval
from sklearn.svm import SVC
import re

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

# For a peptide link map gives three lists:
# - list of amino acids (1-letter) that might also be matches
# - list of probabilities that these aas are the match
# - list of positions relative to the best match that these aas have
def get_neighborhood_list(linkMap):
    matches = re.findall(r'(\w+(?:-\w+)?)\((\d\.\d{2})\)', linkMap)
    aaList = []; valList = []; posList = []; origPosList = np.arange(len(matches)); valSum = 0
    for i,match in enumerate(matches):
        value = float(match[1])
        if len(match[0])==1:
            aaList += [match[0]]
            valList += [value]
            posList += [origPosList[i]]
        elif match[0][-1].isupper():
            aaList += [match[0][-1]]
            valList += [value]
            posList += [origPosList[i]]
        valSum += value
    valList = np.array(valList) / valSum
    posList -= posList[np.argmax(valList)]
    return aaList, valList, posList

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
nofSamplesPerChunk = 10#10000
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
            X[pairkey] = np.zeros(4)[np.newaxis,:]
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        X[pairkey][0,0] += 1 #number of contributing rows
        X[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
        #aaList, valList, posList = get_neighborhood_list(row['PeptideLinkMap1'])
        #print('ye')
        #print("from pdb:", native.residue_type(row['Start1']-4).name1())
        #print(row['PeptideLinkMap1'])
        #print("from pdb:", native.residue_type(aa1Idx-4).name1())
        #print("from LinkedAminoAcid:", row['Linked AminoAcid 1'])
        #print("pos List in fasta:")
        #for pos in posList:
        #    print(native_fasta[aa1Idx + pos])
        #print(aaList)
        #print(valList)
        #print(posList)
        #    Xtest[pairkey][0,2] += 1 #number of neighborhood appearances
        X[pairkey][0,2] += row['match score'] #score
        X[pairkey][0,-1] = res1pos.distance(res2pos) <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for pairkey in X:
    X[pairkey][0,2] /= X[pairkey][0,0]
for contact in contacts:
    if (contact[0]+4, contact[1]+4) in X:
        X[(contact[0]+4, contact[1]+4)][0,3] = True #contact

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
            Xtest[pairkey] = np.zeros(4)[np.newaxis,:]
        res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
        Xtest[pairkey][0,0] += 1 #number of contributing rows
        Xtest[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
        Xtest[pairkey][0,2] += row['match score'] #score
        Xtest[pairkey][0,-1] = res1pos.distance(res2pos) <= 20 #label
    print("Building test set. Finished chunk number %3d."%chunkCount)

for pairkey in Xtest:
    Xtest[pairkey][0,2] /= Xtest[pairkey][0,0]
for contact in contacts:
    if (contact[0]+4, contact[1]+4) in Xtest:
        Xtest[(contact[0]+4, contact[1]+4)][0,3] = True #contact
        
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
