# Parametrize cages which should be there
from metallicious import supramolecular_structure
import parmed as pmd
#
# cage = supramolecular_structure('Co8L16.pdb', {'Co': (2,4)}, topol='Co8L16.top', LJ_type='merz-opc')
# cage.parametrize(out_coord='out_Co8L16.pdb', out_topol='out_Co8L16.top')
# topol = pmd.load_file('out_Co8L16.top')
# assert len(topol.bonds) == 780
#


cage = supramolecular_structure('Fe5L5.pdb', {'Fe': (2,1)}, topol='Fe5L5.top', LJ_type='merz-opc')
cage.parametrize(out_coord='out_Fe5L5.pdb', out_topol='out_Fe5L5.top')
topol = pmd.load_file('out_Fe5L5.top')
assert len(topol.bonds) == 545



print("AAA")
