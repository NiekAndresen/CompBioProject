##Niek Andresen for Computational Biology Project Summer Term 2017

# This script trains and tests the classifier as described in the report.
# It requires the native distances in a file. Uses a .csv file as input which
# should be created by "removing_rows.py".

# GIVE PATHS:
input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
distance_fname = "./1ao6A.distances"
output_fname = "/home/nieck/HSA_data/results/crossval_weightscans"

import numpy as np
import pandas as pd
import crossval as xval
from sklearn.svm import SVC
import re

# For a peptide link map gives three lists:
# - list of amino acids (1-letter) that might also be matches
# - list of probabilities that these aas are the match (normalized)
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
    if abs(idx1-idx2) < 12: #too low sequence separation
        return False
    if idx1 == idx2:
        return False
    return True
            
#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

# accumulate data about all occuring pairs of AAs in X
print("Reading through the data")
X = dict()
chunkCount = 0
chunks = pd.read_csv(input_fname, chunksize=1e5)
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

# chose training set
trainingSetSize = 15000
trainIdx = np.random.choice(np.arange(len(X)), min(trainingSetSize,len(X)))
Xtrain = np.concatenate([ X[x] for x in [list(X.keys())[i] for i in trainIdx] ], axis=0)
nofPositives = int(np.sum(Xtrain[:,-1]))
#Xtrain = np.concatenate([Xtrain[Xtrain[:,-1]==1], Xtrain[Xtrain[:,-1]==0][:nofPositives]], axis=0) #balance training set
print("training set shape:", Xtrain.shape)
print("proportion of positives in training set:", Xtrain[:,-1].mean())

#crossvalidate and train
classifier = SVC()
classifier.fit(Xtrain[:,:-1], Xtrain[:,-1])
#classifier = xval.cv(Xtrain[:,:-1], Xtrain[:,-1], SVC, {'kernel':['linear','rbf']}, nfolds=5, nrepetitions=2, loss_function=xval.zero_one_loss)#xval.false_discovery_rate) #crossvalidation

# put everything together in one list to then predict it all
pairs = list(X.keys())
data = np.concatenate([X[pair] for pair in pairs], axis=0)
exPred = classifier.decision_function(data[:,:-1])
# chose 1400 best winners
nOfWinnersChosen = min(1400, len(pairs))
winners = np.argpartition(exPred, nofWinnersChosen)[-nOfWinnersChosen:]
# list of (idx pair, prediction, label) for each of the best predictions
winners = [(pairs[winner], exPred[winner], X[pairs[winner]][0,-1]) for winner in winners]

correctMatchesCount = 0
matchesCount = 0
for winner in winners:
    matchesCount += 1
    correctMatchesCount += winner[2]

# output
print('classifier kerneltype:', classifier.kernel)
#print('xval loss (0-1):', classifier.cvloss)
print("precision best winners: %f\n"%(float(correctMatchesCount)/matchesCount))
print("len(X): %d"%len(X))
with open(output_fname, 'w') as f:
    f.write("CLASSIFIER WEIGHT SCANS\n")
    f.write("training set shape: %d %d\n"%(Xtrain.shape[0], Xtrain.shape[1]))
    f.write("classifier type: SVM, kernel: %s\n"%(classifier.kernel))
    #f.write("cvloss: %f\n"%classifier.cvloss)
    f.write("\nmatches found: %d\n"%matchesCount)
    f.write("precision: %f\n"%(float(correctMatchesCount)/matchesCount))

