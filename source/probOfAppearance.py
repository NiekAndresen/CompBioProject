
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
##calculates the overall fraction of contacts that exist in the crystal as well (<8A)
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

with open("/home/niek/HSA_data/Run_identifiers_1_2", 'r') as f:
    runs = f.readline()[:-1].split(',')

##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall fraction of highest scoring contacts that exist in the crystal as well (<8A)
def getHighestScoreFoundRatio(dic):
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
        if dic[key][0][1] < 8.:
            foundCorrect += 1
    if not foundSomeMatch:
        return 0.
    return float(foundCorrect) / count
    
##takes dict where for one experiment there is a list of pairs for each scan.
##calculates the overall fraction of lowest scoring contacts that exist in the crystal as well (<8A)
def getLowestScoreFoundRatio(dic):
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
        if dic[key][-1][1] < 8.:
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

#probability of having a random pair of residues in contact in the crystal
count = 0
success = 0
for aa1Idx in range(1,native.total_residue()+1):
    for aa2Idx in range(1, native.total_residue()+1):
        count += 1
        res1pos = native.residue(aa1Idx).nbr_atom_xyz()
        res2pos = native.residue(aa2Idx).nbr_atom_xyz()
        realdist = res1pos.distance(res2pos)
        if realdist < 100.:
            success += 1
print(float(success)/count)

#probability of having at least one pair of residues in contact if two peptides appear as row
#count = 0
#success = 0
#chunks = pd.read_csv("/home/niek/HSA_data/data_experiment_1_2_nodecoys.csv", usecols=columns, chunksize=1e5)
#for chunk in chunks:
#    for i,row in chunk.iterrows():
#        if(i>1000): break
#        foundOne = False
#        indexPairWarValid = False
#        count += 1
#        try:
#            for aa1Idx in range(int(row['Start1'])-4,int(row['Start1'])-4+int(row['LengthPeptide1'])):
#                if foundOne: break
#                for aa2Idx in range(int(row['Start2'])-4,int(row['Start2'])-4+int(row['LengthPeptide2'])):
#                    if aa1Idx < 1 or aa2Idx < 1: #what AAs are those? I don't have them in the .pdb and .fasta
#                        continue
#                    if aa1Idx > native.total_residue() or aa2Idx > native.total_residue():
#                        continue
#                    indexPairWasValid = True
#                    if foundOne: break
#                    res1pos = native.residue(aa1Idx).nbr_atom_xyz()
#                    res2pos = native.residue(aa2Idx).nbr_atom_xyz()
#                    realdist = res1pos.distance(res2pos)
#                    if realdist < 20.:
#                        success += 1
#                        foundOne = True
#        except ValueError:
#            count -= 1
#            continue
#        if not indexPairWasValid:
#            count -= 1
#print(float(success)/count)



