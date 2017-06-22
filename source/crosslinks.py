import numpy as np
import pandas as pd
import pyrosetta as pr

#native_fname = '/home/nieck/HSA_data/1ao6/1ao6A.pdb'
#fasta_fname = '/home/nieck/HSA_data/1ao6/1ao6A_reconstructed.fasta'
#input_fname = "/home/nieck/HSA_data/crosslinks/20PercentFDR_xiFDR0.csv"

native_fname = '/home/niek/HSA_data/1ao6/1ao6A.pdb'
fasta_fname = '/home/niek/HSA_data/1ao6/1ao6A_reconstructed.fasta'
input_fname = "/home/niek/HSA_data/crosslinks/20PercentFDR_xiFDR0.csv"

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

trueCount = 0
crosslinkCount = 0
for aa1Idx,aa2Idx in zip(fromSite, toSite):
    crosslinkCount += 1
    aa1Idx = int(aa1Idx)
    aa2Idx = int(aa2Idx)
    if aa1Idx < 5 or aa2Idx < 5: #these AAs are not in the .pdb
        continue
    if aa1Idx-4 > native.total_residue() or aa2Idx-4 > native.total_residue():
        continue
    res1pos = native.residue(aa1Idx-4).nbr_atom_xyz()
    res2pos = native.residue(aa2Idx-4).nbr_atom_xyz()
    trueCount += res1pos.distance(res2pos) <= 20

print(float(trueCount) / crosslinkCount) #for 20% FDR: ~50% are below 20Angstrom
