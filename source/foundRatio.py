
import pandas as pd
import numpy as np
import pyrosetta as pr

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
header_fname = '/home/niek/HSA_data/header_reduced'
occuring_experiments_fname = "/home/niek/HSA_data/run_identifiers_1_2"
input_fname = "/home/niek/HSA_data/data_experiment_1_2_reduced.csv"
output_fname = "/home/niek/HSA_data/results/ex_1_2_statistics"

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

##takes dict where for one experiment there is a list of pairs for each scan.
##a pair contains distance and rank. pairs are sorted ascending by rank.
##calculates the overall fraction of contacts that exist in the crystal as well (<thresh)
##taking into account the top <number> of ranks
def getFoundRatio(dic, thresh, number):
    foundCorrect = 0
    count = 0
    if len(dic) == 0:
        return 0.
    foundSomeMatch = False
    for key in dic:
        if len(dic[key])==0:
            continue
        foundSomeMatch = True
        n = min(number, len(dic[key]))
        for match in dic[key][:n]:
            count += 1
            if match[0] < thresh:
                foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count

##takes dict where for one experiment there is a list of pairs for each scan.
##a pair contains distance and rank. pairs are sorted ascending by rank.
##calculates the overall fraction of contacts that exist in the crystal as well (<thresh)
def getTotalFoundRatio(dic, thresh):
    foundCorrect = 0
    count = 0
    if len(dic) == 0:
        return 0.
    foundSomeMatch = False
    for key in dic:
        if len(dic[key])==0:
            continue
        foundSomeMatch = True
        for match in dic[key]:
            count += 1
            if match[0] < thresh:
                foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count
    
##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall correlation between score and distance
def getCorrelation(dic):
    scores = np.array([])
    distances = np.array([])
    if len(dic) == 0:
        return 0.
    foundSomeMatch = False
    for key in dic:
        foundSomeMatch = True
        if len(dic[key])==0:
            continue
        for match in dic[key]:
            scores = np.append(scores, [match[0]])
            distances = np.append(distances, [match[1]])
    if not foundSomeMatch:
        return 0
    return np.corrcoef(scores, distances)[0,1]

with open(occuring_experiments_fname, 'r') as f:
    runs = f.readline()[:-1].split(',')

with open(output_fname, 'w') as out:
    out.write("%-35s  %-7s  %-7s  %-7s  %-7s  %-7s  %-7s  %-7s\n"%("experiment","top1","top2","top3","top5","top10","top20","top all"))

#top match scores
threshold = 20.#Angstrom
#contactFoundRatio = {}#for each experiment the percentage of found matches that correspond to a true contact
for ex in runs:
    scans = {}
    chunks = pd.read_csv(input_fname, usecols=columns, chunksize=1e5)
    for chunk in chunks:
        for i,row in chunk.iterrows():
            if row['Run'] != ex:
                continue
            if not row['Scan'] in scans.keys():
                scans[row['Scan']] = []
            try:
                aa1Idx = int(row['ProteinLink1'])
                aa2Idx = int(row['ProteinLink2'])
            except ValueError:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
                continue
            if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                continue
            res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
            res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
            realdist = res1pos.distance(res2pos)
            #scans[row['Scan']] += [(row['match score'], realdist)]
            scans[row['Scan']] += [(realdist, row['MatchRank'])]
    for scan in scans.keys():
        if not all(scans[scan][i][0] <= scans[scan][i+1][0] for i in range(len(scans[scan])-1)):#if not sorted
            scans[scan].sort(key=lambda match: match[0], reverse=False) #ascending by rank
    with open(output_fname, 'a') as out:
        out.write("%-35s  %07.4f  %07.4f  %07.4f  %07.4f  %07.4f  %07.4f  %07.4f\n"\
            %(ex, getFoundRatio(scans, 20., 1), getFoundRatio(scans, 20., 2), getFoundRatio(scans, 20., 3), getFoundRatio(scans, 20., 5), getFoundRatio(scans, 20., 10), getFoundRatio(scans, 20., 20), getTotalFoundRatio(scans, 20.)))
    #contactFoundRatio[ex] = getFoundRatio(scans, 20.)

#print(contactFoundRatio)




