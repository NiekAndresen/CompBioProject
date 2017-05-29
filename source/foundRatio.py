
import pandas as pd
import numpy as np
import pyrosetta as pr

pr.init()
#load native protein
native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
with open('/home/niek/HSA_data/1ao6/1ao6A.fasta', 'r') as f:
    val = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            val += line.rstrip()
native_fasta = '_____' + val
native = pr.pose_from_pdb(native_fname)

with open('/home/niek/HSA_data/header_reduced', 'r') as hr:
    header = hr.readline()[:-1] #removed \n at the end of line
columns = []
for key in header.split(','):
    columns.append(key)
#replace '_' by ' ' again
columns = list(map(lambda string: string.replace('_',' '), columns))
print(columns)

##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall fraction of contacts that exist in the crystal as well (<thresh)
def getFoundRatio(dic, thresh):
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
            if match[1] < thresh:
                foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count

with open("/home/niek/HSA_data/Run_identifiers_1_2", 'r') as f:
    runs = f.readline()[:-1].split(',')

##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall fraction of highest scoring contacts that exist in the crystal as well (<thresh)
def getHighestScoreFoundRatio(dic, thresh):
    foundCorrect = 0
    count = 0
    if len(dic) == 0:
        return 0.
    foundSomeMatch = False
    for key in dic:
        if len(dic[key])==0:
            continue
        foundSomeMatch = True
        count += 1
        if dic[key][0][1] < thresh:
            foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count
    
##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall fraction of lowest scoring contacts that exist in the crystal as well (<thresh)
def getLowestScoreFoundRatio(dic, thresh):
    foundCorrect = 0
    count = 0
    if len(dic) == 0:
        return 0.
    foundSomeMatch = False
    for key in dic:
        if len(dic[key])==0:
            continue
        foundSomeMatch = True
        count += 1
        if dic[key][-1][1] < thresh:
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

#top match scores
contactFoundRatio = {}#for each experiment the percentage of found matches that correspond to a true contact
for ex in runs:
    scans = {}
    chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)
    for chunk in chunks:
        for i,row in chunk.iterrows():
            if row['Run'] != ex:
                continue
            if not row['Scan'] in scans.keys():
                scans[row['Scan']] = []
            try:
                aa1Idx = int(row['ProteinLink1']) #index in the fasta string from above (preceding '_____')
                aa2Idx = int(row['ProteinLink2']) #offset of ~27/28? or 5?
            except ValueError:
                continue
            if aa1Idx < 5 or aa2Idx < 5: #what AAs are those? I don't have them in the .pdb and .fasta
                continue
            if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
                continue
            res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
            res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
            realdist = res1pos.distance(res2pos)
            scans[row['Scan']] += [(row['match score'], realdist)]
    for scan in scans.keys():
        if not all(scans[scan][i][0] >= scans[scan][i+1][0] for i in range(len(scans[scan])-1)):#if not sorted
            scans[scan].sort(key=lambda match: match[0], reverse=True) #descending by score
    contactFoundRatio[ex] = getFoundRatio(scans, 20.)

print(contactFoundRatio)




