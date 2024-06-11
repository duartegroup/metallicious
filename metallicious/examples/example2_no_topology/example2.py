'''
Parametrization of the homoleptic cage using antechamber interface

Equivalent command line prompt:
metallicious -f ru_pd.xyz -vdw_type uff -metal_and_charges Pd 2 Ru 2 -prepare_topol
'''

from metallicious import supramolecular_structure
cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, LJ_type='uff')
cage.prepare_initial_topology()
cage.parametrize(out_coord='out.pdb', out_topol='out.top')

