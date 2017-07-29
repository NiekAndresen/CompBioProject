##Niek Andresen for Computational Biology Project Summer Term 2017

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import crossval as xval
from sklearn.svm import SVC
import re

athome = False #use directory on my local computer
if not athome:
    fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
    input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
    distance_fname = '/home/nieck/HSA_data/1ao6/1ao6A.distances'
    output_fname = "/home/nieck/HSA_data/results/crossval_weightscans"
else:
    fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
    header_fname = '/home/niek/Computational Biology/CompBioProject/headers/header_reduced'
    distance_fname = '/home/niek/HSA_data/1ao6/1ao6A.distances'
    input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
    output_fname = "/home/niek/HSA_data/results/crossval_weightscans"

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

##checks several properties that a pair of indices needs to have in order to
##be taken into account
def valid_idx_pair(idx1, idx2):
    if idx1 < 5 or idx2 < 5: #these AAs are not in the .pdb
        return False
    if idx1-4 > 578 or idx2-4 > 578: #these AAs are not in the .pdb
        return False
    #if abs(idx1-idx2) < 12: #too low sequence separation
    #    return False
    if idx1 == idx2:
        return False
    return True

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
X = dict()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    for i,row in chunk.iterrows():
        if row['MatchRank'] > 5:
            continue
        aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
        aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
        for i,aa1Idx in enumerate(pos1List):
            for j,aa2Idx in enumerate(pos2List):
                if val1List[i] * val2List[j] < .01: #very unlikely
                    continue
                aa1Idx = int(aa1Idx)
                aa2Idx = int(aa2Idx)
                if not valid_idx_pair(aa1Idx, aa2Idx):
                    continue
                if aa1Idx>aa2Idx: #smaller index first
                    continue #continue, because the swapped order will also occur in the loop
                pairkey = (aa1Idx, aa2Idx)
                if not pairkey in X:
                    X[pairkey] = np.zeros(5)[np.newaxis,:]
                    X[pairkey][0,2] = np.inf #for min()
                X[pairkey][0,0] += val1List[i] * val2List[j] #weighted number of contributing entries
                X[pairkey][0,1] += 1/(row['MatchRank']**2) * val1List[i] * val2List[j] #rank score
                X[pairkey][0,2] = min(X[pairkey][0,1], row['MatchRank']) #minimum rank encoutered
                X[pairkey][0,3] += row['match score'] * val1List[i] * val2List[j] #weighted score
                X[pairkey][0,-1] = dist[(aa1Idx-4, aa2Idx-4)] <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

#chose training set
trainingSetSize = 15000
trainIdx = np.random.choice(np.arange(len(X)), min(trainingSetSize,len(X)))
Xtrain = np.concatenate([ X[x] for x in [list(X.keys())[i] for i in trainIdx] ], axis=0)
nofPositives = int(np.sum(Xtrain[:,-1]))
#Xtrain = np.concatenate([Xtrain[Xtrain[:,-1]==1], Xtrain[Xtrain[:,-1]==0][:nofPositives]], axis=0) #balance training set
print("training set shape:", Xtrain.shape)
print("proportion of positives in training set:", Xtrain[:,-1].mean())

#crossvalidate and train
classifier = xval.cv(Xtrain[:,:-1], Xtrain[:,-1], SVC, {'kernel':['rbf']}, nfolds=5, nrepetitions=2, loss_function=xval.zero_one_loss)#xval.false_discovery_rate)

#go through everything. In each scan take the ones with rank <=5 and let the classifier classify on the data stored in X about this pair.
#take the pairs as winners that win in each scan
winners = list()
chunkCount = 0
chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    chunkCount += 1
    experiments = dict()
    for i,row in chunk.iterrows():
        if row['MatchRank'] > 5:
            continue
        if row['Run'] not in experiments:
            experiments[row['Run']] = dict()
        aa1List, val1List, pos1List = get_neighborhood_list(row['PeptideLinkMap1'], row['Start1'])
        aa2List, val2List, pos2List = get_neighborhood_list(row['PeptideLinkMap2'], row['Start2'])
        for i,aa1Idx in enumerate(pos1List):
            for j,aa2Idx in enumerate(pos2List):
                if val1List[i] * val2List[j] < .01: #very unlikely
                    continue
                if not row['Scan'] in experiments[row['Run']]:
                    experiments[row['Run']][row['Scan']] = list()
                aa1Idx = int(aa1Idx)
                aa2Idx = int(aa2Idx)
                if not valid_idx_pair(aa1Idx, aa2Idx):
                    continue
                if aa1Idx>aa2Idx:
                    aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
                pairkey = (aa1Idx, aa2Idx)
                if not pairkey in experiments[row['Run']][row['Scan']]:
                    experiments[row['Run']][row['Scan']] += [pairkey]
    for ex in experiments:
        for scan in experiments[ex]:
            if(len(experiments[ex][scan]) == 0):
                continue
            pairs = experiments[ex][scan]
            Xex = np.concatenate([X[pair] for pair in pairs], axis=0)
            exPred = classifier.decision_function(Xex[:,:-1])
            maxPred = np.argmax(exPred)
            winners += [(pairs[maxPred], exPred[maxPred], X[pairs[maxPred]][0,-1])]
    print("Finished chunk number %3d."%chunkCount)

#chose the 1400 best winners
nOfWinnersChosen = min(1400, len(winners))
preds = [tup[1] for tup in winners]
bestWinners = [winners[i] for i in np.argpartition(preds, nOfWinnersChosen)[-nOfWinnersChosen:]]
correctMatchesCount = 0
matchesCount = 0
for winner in bestWinners:
    matchesCount += 1
    correctMatchesCount += winner[2]
print("scans precision: %f"%(float(correctMatchesCount)/matchesCount))

nofSamplesPerChunk = 500
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
        for i,aa1Idx in enumerate(pos1List):
            for j,aa2Idx in enumerate(pos2List):
                if val1List[i] * val2List[j] < .01: #very unlikely
                    continue
                aa1Idx = int(aa1Idx)
                aa2Idx = int(aa2Idx)
                if not valid_idx_pair(aa1Idx, aa2Idx):
                    continue
                if aa1Idx>aa2Idx:
                    aa1Idx,aa2Idx = aa2Idx,aa1Idx #smaller index first
                pairkey = (aa1Idx, aa2Idx)
                if not pairkey in Xtest:
                    Xtest[pairkey] = np.zeros(5)[np.newaxis,:]
                    Xtest[pairkey][0,2] = np.inf #for min()
                Xtest[pairkey][0,0] += val1List[i] * val2List[j] #weighted number of contributing entries
                Xtest[pairkey][0,1] += 1/(row['MatchRank']**2) * val1List[i] * val2List[j] #rank score
                Xtest[pairkey][0,2] = min(Xtest[pairkey][0,1], row['MatchRank']) #minimum rank encoutered
                Xtest[pairkey][0,3] += row['match score'] * val1List[i] * val2List[j] #weighted score
                Xtest[pairkey][0,-1] = dist[(aa1Idx-4, aa2Idx-4)] <= 20 #label
    print("Finished chunk number %3d."%chunkCount)

Xtest = np.concatenate([Xtest[x] for x in Xtest], axis=0)

print('classifier kerneltype:', classifier.kernel)
print('xval loss (0-1):', classifier.cvloss)
print("test set shape:", Xtest.shape)
testPredictions = classifier.predict(Xtest[:,:-1])
print("proportion of positives in test set:", Xtest[:,-1].mean())
print("test set discoveries:", testPredictions.sum())
testNofPosPredicted = testPredictions.sum()
testNofFalseDiscoveries = testPredictions[Xtest[:,-1]==0].sum()
print("test set FDR:", float(testNofFalseDiscoveries)/testNofPosPredicted)
print("precision best winners: %f\n"%(float(correctMatchesCount)/matchesCount))
print("len(X): %d"%len(X))
with open(output_fname, 'w') as f:
    f.write("CLASSIFIER WEIGHT SCANS\n")
    f.write("training set shape: %d %d\n"%(Xtrain.shape[0], Xtrain.shape[1]))
    f.write("classifier type: SVM, kernel: %s\n"%(classifier.kernel))
    f.write("test set shape: %d %d\n"%(Xtest.shape[0], Xtest.shape[1]))
    f.write("test set discoveries: %d\n"%int(testPredictions.sum()))
    f.write("cvloss: %f\n"%classifier.cvloss)
    f.write("test precision: %f\n"%(1-float(testNofFalseDiscoveries)/testNofPosPredicted))
    f.write("\nmatches found: %d\n"%matchesCount)
    f.write("precision: %f\n"%(float(correctMatchesCount)/matchesCount))

