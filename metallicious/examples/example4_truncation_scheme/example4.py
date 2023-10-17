'''
Parametrization cage using trunction scheme

Equivalent command line prompt:
metallicious -f cage.pdb -p topol.top -vdw_type merz-opc -metal_and_charges Pd 2 -truncate dihedral
'''

import time
start = time.time()
from metallicious import supramolecular_structure
# This will not work becasue there is no exact template for this site:
# cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')
cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc', truncation_scheme='dihedral')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')

print("time", time.time()-start)