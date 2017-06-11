
import pandas as pd
import numpy as np
import pyrosetta as pr
import matplotlib.pyplot as plt

native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/nieck/CompBioProject/headers/header_reduced'
input_fname = "/home/nieck/HSA_data/SDA_HSA_Everything_reduced.csv"
fig_fname = "/home/nieck/HSA_data/results/new_score_precision_everything.png"

#native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
#fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
#header_fname = '/home/niek/HSA_data/header_reduced'
#input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
#fig_fname = "/home/niek/HSA_data/results/new_score_precision_ex_1_2.png"

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

newScore = {} #storing the new score for every pair of linkage sites
averageScore = {} #storing the average score for every pair of linkage sites

chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
for chunk in chunks:
    for i,row in chunk.iterrows():
        site1 = int(row['ProteinLink1'])
        site2 = int(row['ProteinLink2'])
        pairkey = (site1, site2)
        if not pairkey in newScore:
            newScore[pairkey] = 0
            averageScore[pairkey] = [0., 0] #sum of scores and number of occurrences (to calculate average)
        newScore[pairkey] += 1/(row['MatchRank']**2)
        averageScore[pairkey][0] += row['match score']
        averageScore[pairkey][1] += 1
for key in averageScore:
    averageScore[key] = averageScore[key][0] / averageScore[key][1]

def true_matches(scoreDict):
    #sort
    pairList = np.array(list(scoreDict.keys()))
    scoreList = np.array(list(scoreDict.values()))
    sortedIdx = np.argsort(scoreList)
    sortedPairs = pairList[sortedIdx]
    trueMatches = np.zeros(len(sortedPairs))
    for i,pair in enumerate(sortedPairs):
        if pair[0] < 5 or pair[1] < 5: #these AAs are not in the .pdb
            continue
        if pair[0]-4 > native.total_residue() or pair[1]-4 > native.total_residue():
            continue
        res1pos = native.residue(int(pair[0]-4)).nbr_atom_xyz()
        res2pos = native.residue(int(pair[1]-4)).nbr_atom_xyz()
        realdist = res1pos.distance(res2pos)
        trueMatches[i] = realdist < 25 and realdist > 15
    return trueMatches

trueMatches_newScore = true_matches(newScore)
trueMatches_averageScore = true_matches(averageScore)

x_newScore = np.arange(len(trueMatches_newScore))
x_averageScore = np.arange(len(trueMatches_averageScore))
y_newScore = np.cumsum(trueMatches_newScore) / x_newScore
y_averageScore = np.cumsum(trueMatches_averageScore) / x_averageScore
plt.figure()
plt.plot(x_newScore,y_newScore, label="new score")
plt.plot(x_averageScore,y_averageScore, label="average score")
plt.xlabel("number of top scores used")
plt.ylabel("precision")
plt.legend()
plt.savefig(fig_fname)
plt.show()
