import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import crossval as xval
from sklearn.svm import SVC
import re

athome = True #use directory on my local computer
if not athome:
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
    distance_fname = '/home/nieck/HSA_data/1ao6/1ao6A.distances'
    output_fname = "/home/nieck/HSA_data/results/crossval"
else:
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
    distance_fname = '/home/niek/HSA_data/1ao6/1ao6A.distances'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
    output_fname = "/home/niek/HSA_data/results/crossval"

# For a peptide link map gives three lists:
# - list of amino acids (1-letter) that might also be matches
# - list of probabilities that these aas are the match
# - list of index positions in the protein
def get_neighborhood_list(linkMap, startIdx):
    matches = re.findall(r'(\w+(?:-\w+)?)\((\d\.\d{2})\)', linkMap)
    aaList = []; valList = []; posList = []; origPosList = np.arange(len(matches)); valSum = 0
    for i,match in enumerate(matches):
        value = float(match[1])
        if len(match[0])==1:
            aaList += [match[0]]
            valList += [value]
            posList += [origPosList[i] + startIdx]
        elif match[0][-1].isupper():
            aaList += [match[0][-1]]
            valList += [value]
            posList += [origPosList[i] + startIdx]
        valSum += value
    valList = np.array(valList) / valSum
    return aaList, valList, posList

#load fasta
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()
            
#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

#get column names from header
with open(header_fname, 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

#take a number of rows out of each experiment
nofSamplesPerChunk = 5000#10000
X = dict()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    rows = chunk.sample(n=nofSamplesPerChunk)
    for i,row in rows.iterrows():
        if row['MatchRank'] > 5:
            continue
        aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
        aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
        for i, aa1Idx, aa2Idx in zip(range(len(pos1List)), pos1List, pos2List):
            aa1Idx = int(aa1Idx)
            aa2Idx = int(aa2Idx)
            if aa1Idx == aa2Idx:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                continue
            if aa1Idx-4 > 578 or aa2Idx-4 > 578:
                continue
            if aa1Idx>aa2Idx:
                aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
            pairkey = (aa1Idx, aa2Idx)
            if not pairkey in X:
                X[pairkey] = np.zeros(6)[np.newaxis,:]
                X[pairkey][0,2] = np.inf #for min()
            X[pairkey][0,0] += 1 #number of contributing entries
            X[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
            X[pairkey][0,2] = min(X[pairkey][0,1], row['MatchRank']) #minimum rank encoutered
            X[pairkey][0,3] += row['match score'] * val1List[i] * val2List[i] #weighted score
            X[pairkey][0,4] += row['match score'] #score
            X[pairkey][0,-1] = dist[(aa1Idx-4, aa2Idx-4)] <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for pairkey in X: #make score to average score
    X[pairkey][0,4] /= X[pairkey][0,0]

X = np.concatenate([X[x] for x in X], axis=0)
nofPositives = int(np.sum(X[:,-1]))
X = np.concatenate([X[X[:,-1]==1], X[X[:,-1]==0][:nofPositives]], axis=0)
print("training set shape:", X.shape)
print("proportion of positives in training set:", X[:,-1].mean())

#get a test set
Xtest = dict()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    rows = chunk.sample(n=nofSamplesPerChunk)
    for i,row in rows.iterrows():
        if row['MatchRank'] > 5:
            continue
        aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
        aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
        #print("ye", aa1List[i], native.residue_type(int(pos1List[i]-4)).name1())
        for i, aa1Idx, aa2Idx in zip(range(len(pos1List)), pos1List, pos2List):
            aa1Idx = int(aa1Idx)
            aa2Idx = int(aa2Idx)
            if aa1Idx == aa2Idx:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                continue
            if aa1Idx-4 > 578 or aa2Idx-4 > 578:
                continue
            if aa1Idx>aa2Idx:
                aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
            pairkey = (aa1Idx, aa2Idx)
            if not pairkey in Xtest:
                Xtest[pairkey] = np.zeros(6)[np.newaxis,:]
                Xtest[pairkey][0,2] = np.inf #for min()
            Xtest[pairkey][0,0] += 1 #number of contributing entries
            Xtest[pairkey][0,1] += 1/(row['MatchRank']**2) #rank score
            Xtest[pairkey][0,2] = min(Xtest[pairkey][0,1], row['MatchRank']) #minimum rank encoutered
            Xtest[pairkey][0,3] += row['match score'] * val1List[i] * val2List[i] #weighted score
            Xtest[pairkey][0,4] += row['match score'] #score
            Xtest[pairkey][0,-1] = dist[(aa1Idx-4, aa2Idx-4)] <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

for pairkey in Xtest: #make score to average score
    Xtest[pairkey][0,4] /= Xtest[pairkey][0,0]
        
Xtest = np.concatenate([Xtest[x] for x in Xtest], axis=0)

#crossvalidate and train
classifier = xval.cv(X[:,:-1], X[:,-1], SVC, {'kernel':['linear', 'rbf']}, nfolds=5, nrepetitions=2, loss_function=xval.zero_one_loss)#xval.false_discovery_rate)
print('classifier kerneltype:', classifier.kernel)
print('xval loss (0-1):', classifier.cvloss)

print("test set shape:", Xtest.shape)
testPredictions = classifier.predict(Xtest[:,:-1])
print("proportion of positives in test set:", Xtest[:,-1].mean())
print("test set discoveries:", testPredictions.sum())
testNofPosPredicted = (testPredictions==1).sum()
testNofFalseDiscoveries = (np.logical_and(testPredictions==1, Xtest[:,-1]==0)).sum()
print("test set FDR:", float(testNofFalseDiscoveries)/testNofPosPredicted)
with open(output_fname, 'w') as f:
    f.write("training set shape: %d %d\n"%(X.shape[0], X.shape[1]))
    f.write("classifier type: SVM, kernel: %s\n"%(classifier.kernel))
    f.write("test set shape: %d %d\n"%(Xtest.shape[0], Xtest.shape[1]))
    f.write("test set discoveries: %d\n"%int(testPredictions.sum()))
    f.write("cvloss: %f\n"%classifier.cvloss)
