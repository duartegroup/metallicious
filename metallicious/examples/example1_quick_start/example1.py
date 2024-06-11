'''
Parametrization of the cage, using non-bonded topology of the whole cage

Equivalent bash command:
metallicious -f ru_pd.pdb -p ru_pd.top -vdw_type uff -metal_and_charges Pd 2 Ru 2
'''

from metallicious import supramolecular_structure
cage = supramolecular_structure('ru_pd.pdb', metal_charges={'Ru': 2, 'Pd':2 }, topol='ru_pd.top', LJ_type='uff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
