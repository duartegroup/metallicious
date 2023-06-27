import argparse
# TODO this should be class not INFO file

from seminario import single_seminario
from charges import calculate_charges2
from initial_site import create_initial_topol2
from copy_topology_params import copy_bonds, copy_angles, copy_dihedrals

import parmed as pmd
import numpy as np
import MDAnalysis
import os
from load_fingerprint import guess_fingerprint

