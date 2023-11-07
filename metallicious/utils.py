
import os
import shutil
import re

from rdkit import Chem
from tempfile import mkdtemp
import numpy as np

def new_directory(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)

def strip_numbers_from_atom_name(atom_name):
    return re.match("([a-zA-Z]+)", atom_name).group(0)

def strip_numbers_from_atom_names(atom_names):
    return np.array([strip_numbers_from_atom_name(atom_name) for atom_name in atom_names])

