
import metallicious
from metallicious.parametrize_new_sites import supramolecular_structure

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

# print("A")
#
# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Co_Lusby/")
# cage = supramolecular_structure('cage.pdb', {'Co': (3,1)}, vdw_type='merz-tip3p')
# cage.prepare_initial_topology()
# cage.parametrize()

#
# cage = supramolecular_structure(filename='complex.pdb', topol='complex.top', metal_charges={'Pd': 2},
#                                 vdw_type='merz-opc3')
# cage.parametrize(out_coord='out.pdb', out_topol='out.top')
#
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Pd_Lusby/")
# cage = supramolecular_structure('complex.pdb', {'Pd': (2, 1)}, topol='complex.top', vdw_type='merz-opc') # TODO if vdw_type is specified it should overwrite the old one
# cage.prepare_initial_topology()

#cage.parametrize()
#
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Pd_Lusby/")
#
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',  metal_charge_mult={'Pd': (2, 1)},
#                                 topol='init_topol/noncovalent_complex.top', vdw_type='merz-opc', search_library=False)
# cage.parametrize()



#os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/homoleptic/14/temp/")

# os.chdir("/u/fd/chem1540/Downloads/charmm-gui-8906621029/gromacs/")
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',topol='init_topol/noncovalent_complex.top',
#                                 metal_charges={'Pd': 2},  vdw_type='uff', truncation_scheme='dihedral')
# cage.parametrize()

#
# os.chdir("/u/fd/chem1540/temp/temp")
# cage = supramolecular_structure('3_concat.pdb', metal_charge_mult={'Fe': (2,1) },  vdw_type='merz-opc')
# cage.prepare_initial_topology()
# cage.parametrize()
# #
#
#
# os.chdir("/u/fd/chem1540/temp/temp")
# cage = supramolecular_structure('start.pdb', metal_charge_mult={'Fe': (2,1) },  vdw_type='merz-opc')
# cage.prepare_initial_topology(homoleptic_ligand_topol='linker.itp')
# cage.parametrize()


#cage = supramolecular_structure('cage.xyz', {'Pd': (2, 1)})
#cage.prepare_initial_topology()
# cage.parametrize()
#
# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_04/Ga_cage')
# # # cage = supramolecular_structure('cage.xyz', {'Ga': (3,1)}, vdw_type='uff')
# # # cage.prepare_initial_topology(method='gaff', homoleptic_ligand_topol='ligand.top')
# cage = supramolecular_structure('nonbonded.pdb', {'Ga': (3, 1)}, topol = 'noncovalent_complex.top', vdw_type='uff')
# cage.extract_unique_metal_sites()
# cage.unique_sites[0].ligand_charges = [-2, -2, -2]
# cage.parametrize()
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/wierdo")
# cage = supramolecular_structure('wierdo.pdb', {'Pd': (2, 1)}, vdw_type='uff')
# cage.prepare_initial_topology()
# cage.extract_unique_metal_sites()
# print("A")
#
#
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Pd_fujita/")
# # cage = supramolecular_structure('fujita.xyz', {'Pd': (2, 1)}, vdw_type='merz-opc')
# # cage.prepare_initial_topology()
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', {'Pd': (2, 1)}, topol = 'init_topol/noncovalent_complex.top', vdw_type='merz-opc')
# cage.parametrize()

# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Co_cage/")
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', {'Co': (2, 4)},
#                                 topol='init_topol/noncovalent_complex.top', vdw_type='merz-tip3p')
# #cage.prepare_initial_topology()
# cage.extract_unique_metal_sites()
# for site in cage.unique_sites:
#     site.fp_coord_file = f"{site.directory:}/template.pdb"
#     site.fp_topol_file = f"{site.directory:}/template.top"
# cage.assign_fingerprints()
# cage.parametrize()

# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Ru_Pd_cage")
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', {'Ru': (2,1), 'Pd':(2,1)},
#                                 topol='init_topol/noncovalent_complex.top', vdw_type='merz-opc', search_library=False)
# cage.parametrize()
#
# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Zn_knot')
# cage = supramolecular_structure('knot.gro', {'Zn': (2,1)}, vdw_type='merz-opc', improper_metal=False)
# cage.prepare_initial_topology()
# #cage.parametrize()
#
#
# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Fe_cage')
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',{'Fe': (2, 1)},
#                                 topol='init_topol/noncovalent_complex.top', vdw_type='merz-opc', search_library=False)
# cage.parametrize()

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_wb97/Fe_cage')
# #cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',{'Fe': (2, 1)}, topol='init_topol/noncovalent_complex.top', vdw_type='merz-opc', search_library=False)
# cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', {'Fe': (2,1)},
#                                 topol='init_topol/noncovalent_complex.top', vdw_type='merz-opc',
#                                 keywords=['wB97X-D3', 'def2-TZVPP', 'tightOPT', 'freq'], search_library=False)
#
# cage.parametrize()



# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Fe_knot')
# cage = supramolecular_structure('noncovalent_complex.pdb',{'Fe': (2, 1)}, vdw_type='merz-opc')
# cage.prepare_initial_topology()
# #cage.extract_unique_metal_sites()
# #cage.prepare_initial_topology()
# cage.parametrize()

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_06/Fe_knot2')
# cage = supramolecular_structure('noncovalent_complex.pdb', {'Fe': (2, 1)}, vdw_type='merz-opc')
# cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')
# cage.parametrize()
# #
# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Zn_protein")
# cage = supramolecular_structure("protein.gro", topol="protein.top", metal_charge_mult={'Zn': (2, 1)}, improper_metal=False)
# # #cage.prepare_initial_topology()
# cage.parametrize()

os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Co_mof")
name = "site_Co0"
metal_name = "Co"
library_name = "merz-opc"

cage = supramolecular_structure(f'{name:}/saturated_template.xyz', metal_charges={metal_name: 2}, vdw_type=library_name,
                                search_library=False, truncation_scheme='bond', covalent_cutoff=2.5)

cage.prepare_initial_topology()
cage.extract_unique_metal_sites()
cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
cage.sites[0].fp_topol_file = f'{name:}/template.top'
cage.sites[0].load_fingerprint()
cage.sites[0].set_cutoff()
cage.unique_sites = []
cage.parametrize(out_coord='saturated_template_bond.pdb', out_topol='saturated_template_bond.top')

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Zn_mof')
# cage = supramolecular_structure('topol/enormous.gro', {'Zn': (2,1)}, topol='topol/enormous.top', vdw_type='merz-opc')
# cage.parametrize(out_topol='out_enormous.top', out_coord='out_enormous.pdb')

# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Co_mof')
# cage = supramolecular_structure('topol/enormous.gro', {'Co': (2,4)}, topol='topol/enormous.top', vdw_type='merz-opc')
# cage.parametrize(out_topol='out_enormous.top', out_coord='out_enormous.pdb')
#
# #
# os.chdir('/u/fd/chem1540/Research/DA/pd2l4_6+/mix_cn/3_3/MD')
# cage = supramolecular_structure('reopt.xyz', {'Pd': (2,1)}, vdw_type='merz-opc', truncation_scheme='angle')
# #cage.prepare_initial_topology(homoleptic_ligand_topol='linker_top/topol.top') # TODO
# cage.prepare_initial_topology(homoleptic_ligand_topol='new_linker.top')
# cage.parametrize()
# print("A")




# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_06/Co_mof')
# cage = supramolecular_structure('Co_tet.pdb', {'Co': (2,4)}, vdw_type='merz-tip3p',  improper_metal=False)
# cage.prepare_initial_topology(homoleptic_ligand_topol='linker0_new.top')
#cage.parametrize()

'''
os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Cu_protein')
cage = supramolecular_structure('protein.pdb', {'Cu': (2,2), 'Co':(2,1)}, vdw_type='merz-tip3p')
cage.parametrize()
'''

#
# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Ga_cage')
# # cage = supramolecular_structure('cage.xyz', {'Ga': (3,1)}, vdw_type='uff')
# # cage.prepare_initial_topology(method='gaff', homoleptic_ligand_topol='ligand.top')
# cage = supramolecular_structure('nonbonded.pdb', {'Ga': (3, 1)}, topol = 'noncovalent_complex.top', vdw_type='uff')
# cage.extract_unique_metal_sites()
# cage.unique_sites[0].ligand_charges = [-2, -2, -2]
# cage.parametrize()

#0.25105526
#0.0720066
# 0.2357318    0.027520218

#
# cage = supramolecular_structure(filename='noncovalent_complex.pdb', topol='noncovalent_complex.top', metal_charge_mult={'Fe': (2, 1)}, vdw_type='merz-opc3')
# #cage.prepare_initial_topology(method='gaff')
# cage.parametrize(out_coord='out.pdb', out_topol='out.top')


# cage.prepare_initial_topology()
# cage.parametrize()

'''

'''



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

'''

'''
os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Zn_protein/gromacs")
cage = supramolecular_structure('gmx_corrected.pdb', {'Zn': (2, 1)}, vdw_type='merz-tip3p', topol='temp.top')
cage.parametrize()
'''
# #
# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages/Se/cage")
# cage = supramolecular_structure('out.pdb', {'Pd': (2, 1)}, vdw_type='merz-opc', topol='temp.top')
# cage.prepare_initial_topology( homoleptic_ligand_topol='ligand.itp')
# cage.parametrize()

# os.chdir("/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_06_07/Ru_Pd_cage")


# os.chdir("/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_TZVPP/Ru_Pd_cage")
# cage = supramolecular_structure('ru_pd.xyz', {'Ru': (2, 1), 'Pd':(2, 1)}, vdw_type='uff', search_library=False)
# cage.prepare_initial_topology()
# cage.parametrize()


#/u/fd/chem1540/Downloads/charmm-gui-8908552491/gromacs
#os.chdir("/u/fd/chem1540/Downloads/charmm-gui-8908552491/gromacs")

#os.chdir("/u/fd/chem1540/Downloads/charmm-gui-8983783989/gromacs")
# import sys
# home_path = '/u/fd/chem1540'
# sys.path.insert(0,home_path + '/github/cgbind')
# from cgbind import Linker, Cage, ob
# import MDAnalysis
# #syst = MDAnalysis.Universe("step3_input.pdb")
# #syst.atoms.write("linker.xyz")
# linker = Linker(smiles="C1=[N]C=CC(=C1)C2=CC(=CC(=C2)[N+]3=CC=C(C=C3)C4=CC=[N+](C=C4)C)C5=CC=[N]C=C5", name='m4l6_linker', arch_name='m12l24')
# cage = Cage(linker, metal='Pd', metal_charge='2')
# cage.optimise(ob)
# cage.print_xyz_file(filename='cage.xyz')



#
# os.chdir("/u/fd/chem1540/Research/flexiblity/")
# cage = supramolecular_structure(filename='init_topol_weirdo/noncovalent_complex.pdb',
#                                 topol='init_topol_weirdo/noncovalent_complex.top',
#                                 metal_charges={'Pd': 2}, vdw_type='merz-opc', truncation_scheme='angle')
# cage.parametrize()


# cage = supramolecular_structure('cage.xyz', {'Pd': (2, 1)},  vdw_type='merz-tip3p')
# cage.prepare_initial_topology(homoleptic_ligand_topol='LA2.top')
# cage.parametrize()


