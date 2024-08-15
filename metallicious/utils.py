
import os
import re

import MDAnalysis
import rdkit
import numpy as np

def new_directory(directory):
    '''
    Creates new directory
    :param directory: (str) name of new directory
    :return:
    '''
    if not os.path.isdir(directory):
        os.mkdir(directory)
        return True
    else:
        return False


def strip_numbers_from_atom_name(atom_name):
    '''
    Strips the numbers from the string (e.g., C4 -> C)
    :param atom_name: (str) atom name with possible number
    :return: (str)
    '''
    return re.match("([a-zA-Z]+)", atom_name).group(0)

def strip_numbers_from_atom_names(atom_names):
    '''
    Strips numbers from the list of strings, useful for atom names from PDB

    :param atom_names: (list(str))
    :return: (list(str))
    '''
    return np.array([strip_numbers_from_atom_name(atom_name) for atom_name in atom_names])

def guess_aromaticities(syst):
    '''
    Find aromaticity using MDAnalysis, but skip DUMMY atoms as they break the code.
    :param syst: (MDAnalysis.Universe.atoms) input structure
    :return: (list(bool)) list of aromaticy of all atoms
    '''
    aromaticity = [False] * len(syst)
    not_dummy = [idx for idx, name in enumerate(syst.elements) if name != 'DUMMY']

    #aromaticity_from_mdanalysis = MDAnalysis.topology.guessers.guess_aromaticities(syst[not_dummy], force=True)
    mol = syst[not_dummy].convert_to("RDKIT", force=True)
    aromaticity_from_mdanalysis = np.array([atom.GetIsAromatic() for atom in mol.GetAtoms()])

    for idx, _ in enumerate(not_dummy):
        aromaticity[not_dummy[idx]] = aromaticity_from_mdanalysis[idx]
    return aromaticity

def guess_chirality(syst):
    '''
    Returns chirality of the atoms using rdkit

    :param syst: (MDAnalysis.Universe.atoms) input structure
    :return: (list(str)) list of chirality of all atoms
    '''
    chirality = [rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED] * len(syst)
    not_dummy = [idx for idx, name in enumerate(syst.elements) if name != 'DUMMY']
    mol1 = syst[not_dummy].convert_to("RDKIT", force=True)

    for idx, atom in enumerate(mol1.GetAtoms()):
        chirality[not_dummy[idx]]  = atom.GetChiralTag()

    return chirality
