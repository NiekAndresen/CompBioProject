import numpy as np
import pandas as pd
import pyrosetta as pr

#native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
#fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
#input_fname = "/home/nieck/HSA_data/crosslinks/20PercentFDR_xiFDR0.csv"

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
input_fname = "/home/niek/HSA_data/crosslinks/20PercentFDR_xiFDR0.csv"
distance_fname = '/home/niek/HSA_data/1ao6/1ao6A.distances'

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

crosslinks = pd.read_csv(input_fname)
fromSite = crosslinks['fromSite']
toSite = crosslinks['ToSite']
#print("maximum index:", max(fromSite.max(), toSite.max()))
#print("minimum index:", min(fromSite.min(), toSite.min()))

offset = 28
#res1pos = native.residue(186-offset).nbr_atom_xyz()#158
#res2pos = native.residue(159-offset).nbr_atom_xyz()#131
#print(res1pos.distance(res2pos)) #5.878772 is correct

#read native distances
dist = dict()
with open(distance_fname, 'r') as d:
    for line in d:
        arr = line[:-1].split()
        dist[(int(arr[0]), int(arr[1]))] = float(arr[2])

trueCount = 0
crosslinkCount = 0
for aa1Idx,aa2Idx in zip(fromSite, toSite):
    crosslinkCount += 1
    aa1Idx = int(aa1Idx)
    aa2Idx = int(aa2Idx)
    if aa1Idx <= offset or aa2Idx <= offset: #these AAs are not in the .pdb
        continue
    if aa1Idx-offset > native.total_residue() or aa2Idx-offset > native.total_residue():
        continue
    res1pos = native.residue(aa1Idx-offset).nbr_atom_xyz()
    res2pos = native.residue(aa2Idx-offset).nbr_atom_xyz()
    trueCount += res1pos.distance(res2pos) <= 20

print("precision:", float(trueCount) / crosslinkCount) #for 20% FDR: ~67% are below 20Angstrom
