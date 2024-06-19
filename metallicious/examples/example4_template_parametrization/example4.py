'''
Parametrization cage using trunction scheme

Equivalent command line prompt:
metallicious -f cage.pdb -p topol.top -vdw_type merz-opc -metal_and_charges Pd 2 -truncate dihedral
'''

from metallicious import supramolecular_structure
# This will not work becasue there is no exact template for this site:
# cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')

cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges_mult={'Pd': (2,1)}, LJ_type='merz-opc')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')

