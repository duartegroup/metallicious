'''
Parametrization of the cage, using non-bonded topology of the whole cage

Equivalent bash command:
metallicious -f protein.gro -p protein.top -vdw_type merz-opc -metal_and_charges Zn 2
'''

from metallicious import supramolecular_structure
cage = supramolecular_structure('protein.gro', metal_charges={'Zn': 2}, topol='protein.top', LJ_type='merz-opc')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
