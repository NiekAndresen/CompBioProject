
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

chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)

def getFoundRatio(dic):
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
            if match[1] < 8.:
                foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count

#count number of different experiments
runs = set()
for chunk in chunks:
    for i,row in chunk.iterrows():
        runs.add(row['Run'])
nofruns = len(runs)

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
                aa2Idx = int(row['ProteinLink2'])
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
        scans[scan] = sorted(scans[scan], key=lambda match: match[0], reverse=True)#descending by score
    contactFoundRatio[ex] = getFoundRatio(scans)

print(contactFoundRatio)




