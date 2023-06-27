import sys

sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')

from parametrize_new_sites import supramolecular_structure

import os

# os.chdir('/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages/protein')
# os.chdir('/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_05_30/Zn_protein/')
'''
os.chdir('/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_05_30/Ru_Pd_cage/')

cage = supramolecular_structure('ru_pd.xyz', {'Ru': (2,1), 'Pd':(2,1)}, vdw_type='uff')
cage.extract_metal_structures()

site = cage.sites[0]
os.chdir(site.filename)
site.create_initial_topol()
site.seminario()

#site.partial_charge()
#site.cut_extra_atoms_and_save().Graph(MDAnalysis.topology.guesser


print("")
os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_05_30/Lusby_cage")

cage = supramolecular_structure('cage.xyz', {'Pd':(2,1)}, vdw_type='merz-tip3p')
#cage.check_if_paramters_available()
cage.parametrize(out_coord='out.pdb', out_topol='out.top')

'''
'''
os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Cd_protein")
cage = supramolecular_structure('protein.pdb', {'Cd': (2, 1)}, vdw_type='merz-tip3p')
cage.parametrize()
# '''
#
# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Co_cage/")
# cage = supramolecular_structure('Co8L16.xyz', {'Co': (3,1)}, vdw_type='merz-tip3p')
# cage.prepare_initial_topology()
# cage.parametrize()
# print("A")
#
# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Co_Lusby/")
# cage = supramolecular_structure('cage.pdb', {'Co': (3,1)}, vdw_type='merz-tip3p')
# cage.parametrize()

#
# cage = supramolecular_structure(filename='complex.pdb', topol='complex.top', metal_charges={'Pd': 2},
#                                 vdw_type='merz-opc3')
# cage.parametrize(out_coord='out.pdb', out_topol='out.top')
#
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Pd_Lusby/")
# cage = supramolecular_structure('complex.pdb', {'Pd': (2, 1)}, topol='complex.top', vdw_type='merz-opc3')
#

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Co_mof')
# cage = supramolecular_structure('Co_tet.pdb', {'Co': (2,4)}, vdw_type='merz-tip3p',  improper_metal=False)
# cage.parametrize()
'''
os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Cu_protein')
cage = supramolecular_structure('protein.pdb', {'Cu': (2,2), 'Co':(2,1)}, vdw_type='merz-tip3p')
cage.parametrize()
'''
#
os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Fe_cage')
#cage = supramolecular_structure('fe.pdb', {'Fe': (2,1)}, vdw_type='merz-tip3p', improper_metal=False)

cage = supramolecular_structure(filename='noncovalent_complex.pdb', topol='noncovalent_complex.top', metal_charge_mult={'Fe': (2, 1)}, vdw_type='merz-opc3')
#cage.prepare_initial_topology(method='gaff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')


# cage.prepare_initial_topology()
# cage.parametrize()

'''
os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Zn_knot')
cage = supramolecular_structure('knot.gro', {'Zn': (2,1)}, vdw_type='merz-tip3p', improper_metal=False, topol='knot.top')
cage.parametrize()
'''

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Fe_knot')
# cage = supramolecular_structure('noncovalent_complex.pdb',{'Fe': (2, 1)}, vdw_type='merz-tip3p', topol='noncovalent_complex.top')
# #cage.extract_unique_metal_sites()
# #cage.prepare_initial_topology()
# cage.parametrize()

#os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Fe_knot')
# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Ga_cage')
# cage = supramolecular_structure('cage.xyz', {'Ga': (3, 1)}, vdw_type='uff')
# cage.prepare_initial_topology(homoleptic_ligand_topol="ligand.top")
# #cage.parametrize()



# cage = supramolecular_structure('protein.pdb', {'Cd': (2,1)}, vdw_type='merz-tip3p')
# cage = supramolecular_structure('Co_tet.pdb', {'Co': (2,4)}, vdw_type='merz-tip3p')
'''
os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Zn_protein")
cage = supramolecular_structure('1sp1.pdb', {'Zn': (2, 1)}, vdw_type='merz-tip3p')
cage.parametrize()
'''
'''
os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Ru_Pd_cage")
cage = supramolecular_structure('ru_pd.xyz', {'Ru': (2,1), 'Pd':(2,1)}, vdw_type='uff')
cage.parametrize()
'''

'''
os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Zn_protein/gromacs")
cage = supramolecular_structure('gmx_corrected.pdb', {'Zn': (2, 1)}, vdw_type='merz-tip3p', topol='temp.top')
cage.parametrize()
'''
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages/Se/cage")
# cage = supramolecular_structure('out.pdb', {'Pd': (2, 1)}, vdw_type='merz-tip3p', topol='temp.top')
# cage.prepare_initial_topology( homoleptic_ligand_topol='ligand.itp')
# cage.parametrize()
