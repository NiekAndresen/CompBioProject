import random
import re
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
from itertools import *
from os.path import *
import os


"""
A collection of modular function to extract, process and assess
residue-residue contacts and constraints.
"""
__author__ = "Mahmoud Mabrouk"
__email__ = "mahmoud.mabrouk@tu-berlin.de"
__version__ = "0.1"
# Path to the pdbtools (see https://code.google.com/p/pdb-tools/)
PDBTOOLS_PATH = "/scratch/mahmoud/tools/pdb_tools/"

def extract_sequence_from_pdb(pdb_fname):
    """
    Extract sequence from pdb and writes it into a fasta file with the same
    name. Wraps around pdb_seq from pdbTools
    Args:
        pdb_fname: Path to the PDB file
    Returns:
        Path to the fasta file containing the sequence
    """
    script = join(PDBTOOLS_PATH, "pdb_seq.py")
    out_fasta_fname = abspath(pdb_fname)[:-3] + "fasta"
    cmd = "python %s %s > %s" % (script, pdb_fname, out_fasta_fname)
    os.system(cmd)
    return out_fasta_fname



def extract_contacts_from_structure(pdb_fname):
    """
    Extracts the native medium- and long-range
    contacts from a protein structure
    Args:
        pdb_fname: Path of the pdb structure
    Returns:
        List of contacts. A contact is represented as
        a tuple of two residues (integer starting with 1!).
        The list is non redundant, that is if (res1, res2)
        is in the list (res2, res1) will not be there.
    """
    ca, cb = extract_dists_from_pdb(pdb_fname, thresh=10)
    native_cst = get_contacts_from_dists(cb)
    return native_cst


def get_contact_accuracy(contacts_fname,
                         in_contacts=None,
                         native_pdb_fname=None,
                         native_cst=None,
                         thresh=8):
    """
    Returns the accuracy of the contacts in the contacts_fname.
    You could either give as input
    i) a non-redundant list of the native contacts
    ii) the path to the native structure.
    Args:
        contacts_fname: Contact file in the CASP format
        native_pdb_fname: Native pdb
    Returns:
        float(accuracy)
    """
    if native_pdb_fname is None and native_cst is None:
        raise Exception("Please provide either the native structure\
                         or the native contacts")
    elif native_cst is None:
        native_cst = load_native_contacts(native_pdb_fname, thresh=thresh)
    if not in_contacts:
        in_cst = read_casp_cst_file(contacts_fname)
    else:
        in_cst = in_contacts

    n_correct_contacts = 0.0
    for cst in in_cst:
        (res1, res2) = cst
        if (res2, res1) in in_cst:
            raise Exception("Check your inputs. The input contacts have\
                             redundancies in them (res1,res2) AND (res2,res1)\
                             are in inputs!")
        if (res1, res2) in native_cst or (res2, res1) in native_cst:
            n_correct_contacts += 1
            #print("correct %s %s" % (res1, res2))
    if len(in_cst) == 0:
        print("Error: No input contacts")
        return 0.0
    accuracy = n_correct_contacts / float(len(in_cst))
    return accuracy


def load_native_contacts(native_pdb_fname, thresh=8):
    ca, cb = extract_dists_from_pdb(native_pdb_fname, thresh=(thresh + 2))
    native_cst = get_contacts_from_dists(cb, thresh=thresh)
    return native_cst


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


def write_casp_cst_file(contacts, contacts_fname):
    """
    Writes a casp contact file from a set of contacts.
    Args:
        contacts: list(tuple(res1, res2))
        contact_fname: Path to the file to be written
    Returns:
        Name of the written contact file
    """
    f = open(contacts_fname, 'w')
    for (res1, res2) in contacts:
        f.write("%s %s 0 8 0.1337\n" % (res1, res2))
    f.close()
    return contacts_fname


def read_casp_cst_file(contact_fname):
    """
    Reads a contact file in the CASP format.
    Returns:
        list(tuple(res1, res2))
    """
    cst = []
    reg = re.compile(r'(\d+)\s+(\d+)\s+.*')
    with open(contact_fname, 'r') as f:
        for line in f:
            m = re.match(reg, line)
            if m:
                cst.append((int(m.groups()[0]), int(m.groups()[1])))
    return cst


def read_fasta(fasta_file):
    file = open(fasta_file, 'r')
    sequence = ''
    for line in file:
        strline = str(line).strip()
        if strline.count(">") == 0:
            sequence += strline
    file.close()
    return sequence


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


def extract_dists_for_res_pairs(pdb_fname, res_pairs):
    """
    Extract all the native distances for the input residue pairs
    from the native structure pdb file.
    Args:
        pdb_fname: Path to the protein native structure
        res_pairs: list of residue pairs
    Returns:
        dict(tuple(res1, res2)-> float distance in ang between CA
        dict(tuple(res1, res2)-> float distance in ang between CB
    """
    dist_const_ca = {}
    dist_const_cb = {}
    struct_copy_ca = {}
    struct_copy_cb = {}
    structure = _load_pdb(pdb_fname)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    struct_copy_ca[residue.get_id()[1]] = ca.get_coord()
                if residue.has_id("CB"):
                    ca = residue["CB"]
                    struct_copy_cb[residue.get_id()[1]] = ca.get_coord()
                else:  # In the case of Glycine
                    ca = residue["CA"]
                    struct_copy_cb[residue.get_id()[1]] = ca.get_coord()
        break  # We only use the first model
    missing_residues = []
    for res1, res2 in res_pairs:
        res1, res2 = int(res1), int(res2)
        # Residue in native structure
        if res1 in struct_copy_ca and res2 in struct_copy_ca:
            dist_const_ca[res1, res2] = _eucl_dist(
                struct_copy_ca[res1], struct_copy_ca[res2])
            dist_const_cb[res1, res2] = _eucl_dist(
                struct_copy_cb[res1], struct_copy_cb[res2])
        else:
            dist_const_ca[res1, res2] = float('NaN')
            dist_const_cb[res1, res2] = float('NaN')
            if not(res1 in struct_copy_ca) and not(res1 in missing_residues):
                print(res1,
                missing_residues.append(res1))
            if not(res2 in struct_copy_ca) and not(res2 in missing_residues):
                print(res2,
                missing_residues.append(res2))
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
    for res1, res2 in cb_dists.iterkeys():
        if _cst(res1, res2) not in contacts:
            if cb_dists[res1, res2] <= thresh:
                seq_sep = abs(res1 - res2)
                if max_seq_sep >= seq_sep >= min_seq_sep:
                    contacts.append(_cst(res1, res2))
    contacts = sorted(contacts, key=lambda x: x[0] + 0.001 * x[1])
    return contacts


def get_len_pdb(pdb_fname):
    """
    Gets the sequential length of a protein from a pdb
    Args:
        pdb_fname: Path to the pdb file
    Returns:
        A Tuple of ints (s1, s2)
        s1: number of residues in the structure
        s2: highest index of a residue in the structure
        Most of the time these two numbers are the same.
        If not, the user should choose which to use.
    """
    n_res = 0
    highest_index = 0
    structure = _load_pdb(pdb_fname)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                n_res += 1
                highest_index = max(highest_index, residue.get_id()[1])
    return n_res, highest_index


def rand_sample_cst(native_cst, accuracy, nbr_cst, protein_len):
    """
    Randomly sample $nbr_cst contacts with an accuracy $accuracy.
    Args:
        correct_cst: List of all native contacts as tuples (res1, res2)
        accuracy: float defined as true_contacts/nbr_cst
        nbr_cst: Number of contacts to sample
        protein_len: Length of the protein
    Return:
        list((res1, res2)) List of contacts
    """
    n_true_contacts = int(accuracy * nbr_cst)
    n_false_contacts = nbr_cst - n_true_contacts
    if n_true_contacts > len(native_cst):
        raise Exception("Error: Please lower the number of\
                         contacts nbr_cst or the accuracy!")
    rand_native_cst = random.sample(native_cst, n_true_contacts)
    all_cst = [c for c in combinations(range(1, protein_len + 1), 2)]
    med_long_all_cst = filter(lambda x: abs(x[0] - x[1]) >= 12, all_cst)
    wrong_med_long_cst = list(set(med_long_all_cst) - set(native_cst))
    rand_wrong_cst = random.sample(wrong_med_long_cst, n_false_contacts)
    results = rand_native_cst + rand_wrong_cst
    results = sorted(results, key=lambda x: x[0] + 0.001 * x[1])
    return results


def contact_wrongness_stats(contacts, native_dists):
    """
    Args:
        contacts: Contact set. list(tuple(res1, res2)).
        native_dists: All native distance in a protein.
    Returns:
        dict
        'rmsd' -> Root Mean Square Distance of wrong contacts
        'packing' -> Average squential distance between wrong contacts
                     devided by the protein length
        'n10' -> (Number, Percentage) of wrong contacts with native distance > 10A
        'n12.5' -> (Number, Percentage of wrong contacts with native distance > 12.5A
        'n15' -> (Number, Percentage) of wrong contacts with native distance > 15A
        'nmedian'-> (Number, Percentage) of wrong contacts with native distance higher
        than the median native distance for long-range contacts in the protein
    """
    results = {}
    dists, seq_len, ress, diffs = [], [], [], []
    for contact in contacts:
        if native_dists[contact] > 8:
            dists.append(native_dists[contact])
            seq_len.append(abs(contact[0] - contact[1]))
            ress.append(contact[0])
            ress.append(contact[1])
    results["mean_dist"] = np.mean(dists)

    results["seq_len"] = np.mean(seq_len)
    dists = np.array(dists)
    results["n10"] = len(dists[dists > 10.0])
    results["n12.5"] = len(dists[dists > 12.5])
    results["n15"] = len(dists[dists > 15.0])
    ress = sorted(ress)
    for i, res in enumerate(ress[:-1]):
        diffs.append((ress[i+1] - res) ** 2)
    results["coverage"] = np.mean(diffs)
    return results


def contact_correctness_stats(contacts, native_dists):
    """
    Args:
        contacts: Contact set. list(tuple(res1, res2))
        native_dists: The CB native distances in the protein
    Returns:
        dict containing stats:
        'mean_dist': Mean distance of correct contacts in A
        'seq_len' : Mean sequence length of correct lengths
        'diff_to_the_mid' : NOT IMPLEMENTED
                            mean distance of correct contacts
                            devided by the median distance between
                            residues in the protein
        'packing' : NOT IMPLEMENTED
                    Average sequential distance between correct
                    contacts devided by the protein length
        'coverage': a stupid measure that shows how much the residues
                    making up the correct contacts differ between each
                    other
    """
    results = {}
    dists, seq_len, ress, diffs = [], [], [], []
    for contact in contacts:
        if native_dists[contact] <= 8:
            dists.append(native_dists[contact])
            seq_len.append(abs(contact[0] - contact[1]))
            ress.append(contact[0])
            ress.append(contact[1])
    results["mean_dist"] = np.mean(dists)

    results["seq_len"] = np.mean(seq_len)

    ress = sorted(ress)
    for i, res in enumerate(ress[:-1]):
        diffs.append((ress[i+1] - res) ** 2)
    results["coverage"] = np.mean(diffs)
    return results


def calc_fuzzy_accuracy(contacts, native_dists, fuzziness=1):
    """
    We define here a contact (res1, res2) as correct if
    there is at least one contact between res1-n..res1+n
    res2-n..res2+n, where n is the fuzziness
    Args:
        contacts: Contact set. list(tuple(res1, res2))
        native_dists: The CB native distances in the protein
    Returns:
        fuzzy accuracy
    """
    acc = 0
    for (res1, res2) in contacts:
        accurate = False
        for fres1 in range(res1-fuzziness, res1+fuzziness+1):
            for fres2 in range(res2-fuzziness, res2+fuzziness+1):
                if (fres1, fres2) in native_dists:
                    if native_dists[(fres1, fres2)] <= 8:
                        accurate = True
        if accurate:
            acc += 1
    return float(acc) / len(contacts)


def compute_contact_satifaction(casp_contacts_fname,
                                native_contacts,
                                ctype="medium+long"):
    """
    Compute fraction of contacts satisfied in a decoy.
    Args:
        casp_contacts_fname: File containing the contacts satisfied in a
                             decoy in the casp format file
        native_contacts: list containing the native contacts
        type: {long, medium, medium+long}
    Returns:
        (float) Fraction of contacts that are satisfied
    """
    contacts = read_casp_cst_file(casp_contacts_fname)
    sat, total = 0, 0
    #print(contacts)
    for c in native_contacts:
        if (ctype == "medium" and
           abs(c[1] - c[0]) >= 12 and
           abs(c[1] - c[0]) < 24):
            total += 1
            if c in contacts:
                sat += 1
        if (ctype == "long" and
           abs(c[1] - c[0]) >= 24):
            total += 1
            if c in contacts:
                sat += 1
        if (ctype == "medium+long" and
           abs(c[1] - c[0]) >= 12):
            total += 1
            if c in contacts:
                sat += 1
    if total > 0:
        frac = float(sat) / total
    else:
        frac = float(-1)
    return frac


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
