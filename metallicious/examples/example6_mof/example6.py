'''
Parametrization of the cage, using non-bonded topology of the whole cage

Equivalent bash command:
metallicious -f protein.gro -p protein.top -vdw_type merz-opc -metal_and_charges Zn 2
'''

from metallicious import supramolecular_structure
cage = supramolecular_structure('Co_mof.pdb', metal_charges={'Co':2}, LJ_type='merz-opc')
cage.prepare_initial_topology()
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
