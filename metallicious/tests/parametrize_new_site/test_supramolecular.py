# Parametrize cages which should be there
from metallicious import supramolecular_structure
import parmed as pmd

def are_bonds_the_same(bond1, bond2):
    if bond1.type == bond2.type and bond1.atom1.idx == bond2.atom1.idx and bond1.atom2.idx == bond2.atom2.idx:
        return True
    else:
        return False
def compare_bonds(topol1, topol2):
    count = 0
    for bond1 in topol1.bonds:
        for bond2 in topol2.bonds:
            if are_bonds_the_same(bond1, bond2):
                count +=1
    return count

def test_patcher():
    cage = supramolecular_structure('Co8L16.pdb', {'Co': (2,4)}, topol='Co8L16.top', LJ_type='merz-opc')
    cage.parametrize(out_coord='out_Co8L16.pdb', out_topol='out_Co8L16.top')
    topol = pmd.load_file('out_Co8L16.top')
    topol_old = pmd.load_file('Co8L16.top')
    assert len(topol.bonds) - len(topol_old.bonds) == 8*6

    cage = supramolecular_structure('Fe5L5.pdb', {'Fe': (2,1)}, topol='Fe5L5.top', LJ_type='merz-opc')
    cage.parametrize(out_coord='out_Fe5L5.pdb', out_topol='out_Fe5L5.top')
    topol = pmd.load_file('out_Fe5L5.top')
    topol_old = pmd.load_file('Fe5L5.top')
    assert len(topol.bonds) - len(topol_old.bonds) == 5*6

    cage = supramolecular_structure('Pd2L4.pdb',  metal_charges={'Pd': 2}, topol='Pd2L4.top', LJ_type='merz-opc')
    cage.parametrize(out_coord='out_Pd2L4.pdb', out_topol='out_Pd2L4.top')
    topol = pmd.load_file('out_Pd2L4.top')
    topol_old = pmd.load_file('Pd2L4.top')
    assert len(topol.bonds) - len(topol_old.bonds) == 2*4

    cage = supramolecular_structure('Pd6Ru8L28.pdb',  metal_charges={'Pd': 2, 'Ru': 2}, topol='Pd6Ru8L28.top', LJ_type='uff')
    cage.parametrize(out_coord='out_Pd6Ru8L28.pdb', out_topol='out_Pd6Ru8L28.top')
    topol = pmd.load_file('out_Pd6Ru8L28.top')
    topol_old = pmd.load_file('Pd6Ru8L28.top')
    assert len(topol.bonds) - len(topol_old.bonds) == 6*4 + 8*6

    cage = supramolecular_structure('Pd6L4.pdb',  metal_charges={'Pd': 2}, topol='Pd6L4.top', LJ_type='merz-opc')
    cage.parametrize(out_coord='out_Pd6L4.pdb', out_topol='out_Pd6L4.top')
    topol = pmd.load_file('out_Pd6L4.top')
    topol_old = pmd.load_file('Pd6L4.top')
    assert len(topol.bonds) - len(topol_old.bonds) == 6*4


def test_truncation_scheme():
    cage = supramolecular_structure('Pd2L4.pdb',  metal_charges={'Pd': 2}, topol='Pd2L4.top',
                                    LJ_type='merz-opc', truncation_scheme='dih')
    cage.parametrize(out_coord='out_Pd2L4.pdb', out_topol='out_Pd2L4.top')
    topol = pmd.load_file('out_Pd2L4.top')
    topol_old = pmd.load_file('Pd2L4.top')
    assert len(topol.bonds) - compare_bonds(topol, topol_old) == 7*4*2

    cage = supramolecular_structure('Pd2L4.pdb',  metal_charges={'Pd': 2}, topol='Pd2L4.top',
                                    LJ_type='merz-opc', truncation_scheme='ang')
    cage.parametrize(out_coord='out_Pd2L4.pdb', out_topol='out_Pd2L4.top')
    topol = pmd.load_file('out_Pd2L4.top')
    topol_old = pmd.load_file('Pd2L4.top')
    assert len(topol.bonds) - compare_bonds(topol, topol_old) == 3*4*2

    cage = supramolecular_structure('Pd2L4.pdb',  metal_charges={'Pd': 2}, topol='Pd2L4.top',
                                    LJ_type='merz-opc', truncation_scheme='bond')
    cage.parametrize(out_coord='out_Pd2L4.pdb', out_topol='out_Pd2L4.top')
    topol = pmd.load_file('out_Pd2L4.top')
    topol_old = pmd.load_file('Pd2L4.top')
    assert len(topol.bonds) - compare_bonds(topol, topol_old) == 1*4*2
