import random
import re
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
from itertools import *
from os.path import *
import os

def extract_dists_from_pdb(pdb_fname, thresh=8):
    """
    Extracts all residue-residue distances (between ca and cb)
    from a pdb structure. Only extract distance lower than thresh.
    In the case of glycine the ca atom will always be used.
    Args:
        pdb_fname: path to the pdb file
        thresh: maximum distance between two residue to be included
                in the results
    Returns:
        dict((res1, res2)->distance in angstrom) (for ca)
        dict((res1, res2)->distance in angstrom) (for cb)
    """
    dist_const_ca = {}
    dist_const_cb = {}
    struct_copy_ca = {}
    struct_copy_cb = {}
    all_residues = []
    structure = _load_pdb(pdb_fname)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                all_residues.append(residue.get_id()[1])
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    struct_copy_ca[residue.get_id()[1]] = ca.get_coord()
                if residue.has_id("CB"):
                    cb = residue["CB"]
                else:  # In the case of Glycine
                    if "CA" in residue:
                        cb = residue["CA"]
                    else:  # Probably a HETATOM
                        continue
                struct_copy_cb[residue.get_id()[1]] = cb.get_coord()
        break  # We only use the first model
    for (res1, res2) in combinations(all_residues, 2):
        res1, res2 = int(res1), int(res2)
        try:
            dist = _eucl_dist(struct_copy_ca[res1], struct_copy_ca[res2])
        except:
            print("Skipping..")
            return {}, {}
        if dist <= thresh:
            dist_const_ca[res1, res2] = _eucl_dist(
                struct_copy_ca[res1], struct_copy_ca[res2])
            dist_const_ca[res2, res1] = dist_const_ca[res1, res2]
            dist_const_cb[res1, res2] = _eucl_dist(
                struct_copy_cb[res1], struct_copy_cb[res2])
            dist_const_cb[res2, res1] = dist_const_cb[res1, res2]
    return dist_const_ca, dist_const_cb

def get_contacts_from_dists(cb_dists, thresh=8, min_seq_sep=12,
                            max_seq_sep=np.Inf):
    """
    Returns a list of residue-residue contacts from
    a dictionary of distances between cb atoms
    Args:
        cb_dists: dict(res1,res2)->float(distance in angstrom)
        thresh: Contact distance threshold
        min_seq_sep: minimal sequence seperation of two residues
                     to be taken into account
        max_seq_sep: maximum sequence seperation of two residues
                     to be taken into account
    Returns:
        list((res1, res2)). Non redundant (if (res1, res2) in list
        then (res2, res1) not in list). res1 < res2.
    """
    contacts = []
    for res1, res2 in cb_dists.keys():
        if _cst(res1, res2) not in contacts:
            if cb_dists[res1, res2] <= thresh:
                seq_sep = abs(res1 - res2)
                if max_seq_sep >= seq_sep >= min_seq_sep:
                    contacts.append(_cst(res1, res2))
    contacts = sorted(contacts, key=lambda x: x[0] + 0.001 * x[1])
    return contacts

def write_rosetta_cst_restraints(contacts, restraint_fname, protein_seq_fname,
                                 contact_type='CB', low_bound=0.0,
                                 upp_bound=8.0, welldepth=2.0):
    """
    Writes a rosetta restraints file from a set of contacts.
    Args:
        contacts: list(tuple(res1, res2))
        restraint_fname: Path to the file to be written
        protein_seq_fname: File containing the protein sequence (to determine
                     Glycine residues)
        contact_type: Between CB or CA atoms
        low_bound: lower bound (0 per default)
        upp_bound: upper bound (8 per default)
        welldepth: Depth for the energy term
    Returns:
        Name of the written restraint file
    """
    protein_seq = read_fasta(protein_seq_fname)
    f = open(restraint_fname, 'w')
    for (res1, res2) in contacts:
        atom_1 = "CA" if _is_glycine(res1, protein_seq) else contact_type
        atom_2 = "CA" if _is_glycine(res2, protein_seq) else contact_type
        f.write("AtomPair %s  %s %s  %s BOUNDED  %s %s %s NOE\n" % (atom_1,
                                                                    res1,
                                                                    atom_2,
                                                                    res2,
                                                                    low_bound,
                                                                    upp_bound,
                                                                    welldepth))
    f.close()
    return restraint_fname

def load_native_contacts(native_pdb_fname, thresh=8):
    ca, cb = extract_dists_from_pdb(native_pdb_fname, thresh=(thresh + 2))
    native_cst = get_contacts_from_dists(cb, thresh=thresh)
    return native_cst
""""""""""""""""""""""""""""""""
""""""" Help Functions """""""""
""""""""""""""""""""""""""""""""


def _is_glycine(residue_n, protein_seq):
    return protein_seq[residue_n - 1] == "G"


def _eucl_dist(point1, point2):
    """Calculate Euclidean distance between two points in 3D space"""
    return sqrt((point1[0] - point2[0]) * (point1[0] - point2[0]) +
                (point1[1] - point2[1]) * (point1[1] - point2[1]) +
                (point1[2] - point2[2]) * (point1[2] - point2[2]))


def _load_pdb(pdbFile):
    p = PDBParser(PERMISSIVE=1)
    structure_id = "native"
    return p.get_structure(structure_id, pdbFile)


def _cst(res1, res2):
    return (min(res1, res2), max(res1, res2))

def read_fasta(fasta_file):
    file = open(fasta_file, 'r')
    sequence = ''
    for line in file:
        strline = str(line).strip()
        if strline.count(">") == 0:
            sequence += strline
    file.close()
    return sequence

""""""""""""""""""""""""""""""""""""
""""""""""""""" main """""""""""""""
""""""""""""""""""""""""""""""""""""
nat_cst = load_native_contacts("/home/niek/HSA_data/1ao6/1ao6A.pdb")
write_rosetta_cst_restraints(nat_cst, "/home/niek/HSA_data/1ao6/1ao6A.cst", "/home/niek/HSA_data/1ao6/1ao6A.fasta")
