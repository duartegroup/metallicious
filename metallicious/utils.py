
import os
import shutil
import re

from rdkit import Chem
from tempfile import mkdtemp
import numpy as np

def new_directory(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)

def mdanalysis_to_rdkit(syst): # TODO delete
    here = os.getcwd()
    tmpdir_path = mkdtemp()
    os.chdir(tmpdir_path)
    syst.atoms.write("mol1.pdb")
    mol = Chem.MolFromPDBFile("mol1.pdb", removeHs=False)
    os.chdir(here)
    shutil.rmtree(tmpdir_path)
    return mol


def strip_numbers_from_atom_names(atom_names):
    return np.array([re.match("([a-zA-Z]+)", atom_name).group(0) for atom_name in atom_names])