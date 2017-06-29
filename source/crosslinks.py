import numpy as np
import pandas as pd

fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
input_fname = "/home/niek/HSA_data/crosslinks/20PercentFDR_xiFDR0.csv"
distance_fname = '/home/niek/HSA_data/1ao6/1ao6A.distances'

#load native protein
with open(fasta_fname, 'r') as f:
    native_fasta = ''
    for line in f:
        if line.startswith('>'):
            continue
        else:
            native_fasta += line.rstrip()

crosslinks = pd.read_csv(input_fname)
fromSite = crosslinks['fromSite']
toSite = crosslinks['ToSite']
fromSite = fromSite[crosslinks['isTT']] #only target-target matches
toSite = toSite[crosslinks['isTT']] #only target-target matches

offset = 28

#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

trueCount = 0
crosslinkCount = 0
for aa1Idx,aa2Idx in zip(fromSite, toSite):
    aa1Idx = int(aa1Idx)
    aa2Idx = int(aa2Idx)
    if aa1Idx <= offset or aa2Idx <= offset: #these AAs are not in the .pdb
        continue
    if aa1Idx-offset > 578 or aa2Idx-offset > 578:
        continue
    if abs(aa1Idx-aa2Idx) < 12: #only higher sequence separation
        continue
    if aa1Idx>aa2Idx:
        aa1Idx,aa2Idx = aa2Idx,aa1Idx
    crosslinkCount += 1
    trueCount += dist[(aa1Idx-offset, aa2Idx-offset)] <= 20
    #result: the two distance measures are the same

print("precision:", float(trueCount) / crosslinkCount) #for 20% FDR: ~60% are below 20 Angstrom, ~71% are below 25 Angstrom
