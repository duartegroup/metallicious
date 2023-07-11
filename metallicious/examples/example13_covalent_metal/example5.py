import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')

import os



#from extract_metal_site import extract_metal_structure
#from seminario import multi_seminario
#from charges import calculate_charges
#from copy_topology_params import copy_bonds, copy_angles, copy_dihedrals
try:
    from metallicious.main import patcher
    from metallicious.extract_metal_site import extract_metal_structure
    from metallicious.initial_site import create_initial_topol

except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure
    from cgbind2pmd.initial_site import create_initial_topol

f='ligand_old.pdb'
metal='Se'
metal_charge=0
vdw_type='uff'

extract_metal_structure(f, metal, "ligand_fp")
topols = create_initial_topol(metal, metal_charge, vdw_data_name=vdw_type)

