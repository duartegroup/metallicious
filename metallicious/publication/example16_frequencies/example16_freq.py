import shutil

from metallicious.parametrize_new_sites import supramolecular_structure
import os

folder_names = ['Ru_Pd_cage', 'Pd_Lusby', 'Co_cage',  'Co_mof',  'Fe_cage',  'Fe_knot',   'Ga_cage',  'Pd_fujita',    'Pd_o', 'Zn_knot',  'Zn_mof']
#folder_names = ['Ga_cage']
uff_folders = ['Ga_cage', 'Ru_Pd_cage' ]

metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
                     'Pd_fujita':2, 'Pd_o':2, 'Ru_Pd_cage':2,   'Zn_knot':2,  'Zn_mof':2,  'Zn_protein':2, 'megacage':2}

#
# for basis_folder in ['10', 'TZVPP', 'wb97']:
#     for idx, _ in list(enumerate(folder_names)):
#         if not (basis_folder=='TZVPP' and folder_names[idx]=='Ru_Pd_cage'):
#             folder_name = folder_names[idx]
#             #os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                             metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                             vdw_type = library_name, search_library=False)
#
#             if folder_name!="Ga_cage":
#                 cage.prepare_initial_topology()
#             else:
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small.top')
#
#             cage.extract_unique_metal_sites()
#             cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#             cage.sites[0].fp_topol_file = f'{name:}/template.top'
#             cage.sites[0].load_fingerprint()
#             cage.sites[0].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord='saturated_template.pdb', out_topol='saturated_template.top')
# folder_names = ['Zn_mof']

# for basis_folder in ['10', 'TZVPP', 'wb97']:
#     for idx, _ in list(enumerate(folder_names)):
#         if not (basis_folder=='TZVPP' and folder_names[idx]=='Ru_Pd_cage'):
#
#             folder_name = folder_names[idx]
#
#             #os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             if folder_name in ['Fe_knot', 'Pd_fujita']:
#                 cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                                 topol='init_topol/noncovalent_complex_resp.top',
#                                                 metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                                 vdw_type=library_name, search_library=False)
#             else:
#                 cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                                 metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                                 vdw_type=library_name, search_library=False)
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small_resp.top')
#
#             cage.extract_unique_metal_sites()
#             cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#             cage.sites[0].fp_topol_file = f'{name:}/template.top'
#             cage.sites[0].load_fingerprint()
#             cage.sites[0].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord='saturated_template_resp.pdb', out_topol='saturated_template_resp.top')

#
# for basis_folder in ['zhang', 'uff']:
#     for idx, _ in list(enumerate(folder_names)):
#         if not (("Pd" in basis_folder or "Ga" in basis_folder) and basis_folder=='zhang'):
#             folder_name = folder_names[idx]
#             #os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10_{basis_folder:}/{folder_name}/")
#
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                             metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                             vdw_type = library_name, search_library=False)
#
#             if folder_name!="Ga_cage":
#                 cage.prepare_initial_topology()
#             else:
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small.top')
#             cage.extract_unique_metal_sites()
#
#             cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#             cage.sites[0].fp_topol_file = f'{name:}/template.top'
#             cage.sites[0].load_fingerprint()
#             cage.sites[0].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord='saturated_template.pdb', out_topol='saturated_template.top')

# basis_folder = '10'
# #
# # folder_names = ['Zn_protein']
# #folder_names = ['Fe_knot']
# folder_names = ['Zn_knot']
#
# # #for basis_folder in ['10']:
# #for truncation_scheme_string in ['dihedral']:
# for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
#     #for truncation_scheme_string in [ 'bond']:
#
#     for idx, folder_name in enumerate(folder_names):
#         if not (folder_name == 'Zn_knot' and (truncation_scheme_string=='dihedral' or truncation_scheme_string=='angle')):
#             print(idx, folder_name, truncation_scheme_string)
#
#             # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                             metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                             vdw_type=library_name, search_library=False, truncation_scheme=truncation_scheme_string, covalent_cutoff=2.5)
#
#             if folder_name != "Ga_cage":
#                 cage.prepare_initial_topology()
#             else:
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small.top')
#             cage.extract_unique_metal_sites()
#
#             cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#             cage.sites[0].fp_topol_file = f'{name:}/template.top'
#             cage.sites[0].load_fingerprint()
#             cage.sites[0].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord=f'saturated_template_{truncation_scheme_string:}.pdb', out_topol=f'saturated_template_{truncation_scheme_string:}.top')
#

#for basis_folder in ['10']:
#for truncation_scheme_string in ['dihedral']:
#folder_names = ['Zn_knot',  'Zn_mof']
#
# basis_folder = '10'
# for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
#     #for truncation_scheme_string in [ 'bond']:
#     for idx, folder_name in enumerate(folder_names):
#         if not (folder_name == 'Zn_knot' and truncation_scheme_string=='angle'):
#             print(idx,folder_name, truncation_scheme_string)
#             # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             if folder_name in ['Fe_knot', 'Pd_fujita']:
#                 cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                                 topol='init_topol/noncovalent_complex_resp.top',
#                                                 metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                                 vdw_type=library_name, search_library=False, truncation_scheme=truncation_scheme_string, covalent_cutoff=2.5)
#             else:
#                 cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                                 metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                                 vdw_type=library_name, search_library=False, truncation_scheme=truncation_scheme_string, covalent_cutoff=2.5)
#
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small_resp.top')
#
#
#             cage.extract_unique_metal_sites()
#
#             cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#             cage.sites[0].fp_topol_file = f'{name:}/template.top'
#             cage.sites[0].load_fingerprint()
#             cage.sites[0].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord=f'saturated_template_{truncation_scheme_string:}_resp.pdb', out_topol=f'saturated_template_{truncation_scheme_string:}_resp.top')


# for idx, _ in list(enumerate(folder_names)):
#     folder_name = folder_names[idx]
#     print(folder_name)
#     os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10_torsion/{folder_name}/")
#
#     name = f"site_{folder_name[:2].title():}0"
#     metal_name = folder_name[:2].title()
#
#     library_name = 'merz-opc'
#     if folder_name in uff_folders:
#         library_name = 'uff'
#
#     cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                     metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                     vdw_type = library_name, search_library=False)
#
#     if folder_name != "Ga_cage":
#         cage.prepare_initial_topology()
#     else:
#         cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small.top')
#     cage.extract_unique_metal_sites()
#
#     cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#     cage.sites[0].fp_topol_file = f'{name:}/template.top'
#     cage.sites[0].load_fingerprint()
#     cage.sites[0].set_cutoff()
#     cage.unique_sites = []
#     cage.parametrize(out_coord='saturated_template_impropers.pdb', out_topol='saturated_template_impropers.top')
#
#
#     cage = supramolecular_structure(f'{name:}/bonded/saturated_template_optimised.xyz',
#                                     metal_charges={metal_name: metal_charges_dic[folder_name]},
#                                     vdw_type = library_name, search_library=False)
#
#     if folder_name != "Ga_cage":
#         cage.prepare_initial_topology()
#     else:
#         cage.prepare_initial_topology(homoleptic_ligand_topol='linker_small.top')
#     cage.extract_unique_metal_sites()
#
#     cage.extract_unique_metal_sites()
#     cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
#     cage.sites[0].fp_topol_file = f'{name:}/template_saturated.top'
#     cage.sites[0].load_fingerprint()
#     cage.sites[0].set_cutoff()
#     cage.unique_sites = []
#     cage.parametrize(out_coord='saturated_template_noimpropers.pdb', out_topol='saturated_template_noimpropers.top')
#

#cage.parametrize(out_coord='saturated_template_bond.pdb', out_topol='saturated_template_bond.top')
#
# folder_names = ['Co_mof',  'Ga_cage', 'Zn_mof',  'Zn_protein']
#
# folder_names_cages = ['Pd_Lusby', 'Co_cage', 'Fe_cage',  'Fe_knot', 'Ga_cage',  'Pd_fujita', 'Pd_o', 'Ru_Pd_cage', 'Zn_knot']
#
# folder_names_noncages = ['Co_mof', 'Zn_mof', 'Zn_protein']
#
# # Pd-O jest specjane bo to jest tak naprese Se, trzeba zrobic pliki topologiczne dla niego
folder_names = ['Pd_Lusby','Co_cage', 'Fe_cage',  'Fe_knot',   'Ga_cage',  'Pd_fujita', 'Zn_knot',  'Ru_Pd_cage', ]
uff_folders = ['Ga_cage', 'Ru_Pd_cage' ]

# metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
#                      'Pd_fujita':2, 'Ru_Pd_cage':2,   'Zn_knot':2,  'Zn_mof':2,  'Zn_protein':2}
#
# metal_improper = {'Pd_Lusby':True,'Co_cage':False, 'Co_mof':False,  'Fe_cage':False,  'Fe_knot':False, 'Ga_cage':False,
#                      'Pd_fujita':True, 'Ru_Pd_cage':True,   'Zn_knot':False,  'Zn_mof':False,  'Zn_protein':False}
#
# #for basis_folder in ['10', 'TZVPP', 'wb97']:

#
# folder_names = ['Co_cage']
#for basis_folder in ['10', 'TZVPP', 'wb97']:
#for basis_folder in ['10_zhang', '10_uff']:
#
# folder_names = ['megacage']
#
#for basis_folder in ['']:
# for basis_folder in ['TZVPP', 'wb97']:
#     for idx, _ in list(enumerate(folder_names)):
#         if not (basis_folder == 'TZVPP' and folder_names[idx] == 'Ru_Pd_cage'):
#             folder_name = folder_names[idx]
#             print(folder_name)
#
#             os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#             if folder_name == 'megacage':
#                 name = f"site_Pd0"
#                 metal_name = 'Pd'
#             else:
#                 name = f"site_{folder_name[:2].title():}0"
#                 metal_name = folder_name[:2].title()
#
#             metal_charges = {metal_name: metal_charges_dic[folder_name]}
#             if folder_name == 'Ru_Pd_cage':
#                 metal_charges = {'Ru' : 2, 'Pd' : 2}
#
#             library_name = 'merz-opc'
#             if folder_name in uff_folders:
#                 library_name = 'uff'
#
#             #shutil.copyfile(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/{folder_name:}/crystal.pdb", "crystal.pdb")
#
#             if folder_name=='megacage':
#                 cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',
#                                                 topol='init_topol/noncovalent_complex.top',
#                                                 metal_charges=metal_charges,
#                                                 vdw_type = library_name, search_library=False)
#             else:
#                 cage = supramolecular_structure(f'crystal.pdb',
#                                                 metal_charges=metal_charges,
#                                                 vdw_type = library_name, search_library=False)
#
#             if folder_name == 'Pd_o' or folder_name == 'megacage':
#                 # shutil.copyfile(
#                 #     f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/{folder_name:}/ligand.itp",
#                 #     "ligand.itp")
#                 cage.prepare_initial_topology(homoleptic_ligand_topol = 'ligand.itp')
#             elif folder_name == 'Ga_cage':
#                 # shutil.copyfile(
#                 #     f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/{folder_name:}/ligand.top",
#                 #     "ligand.top")
#                 cage.prepare_initial_topology(homoleptic_ligand_topol='ligand.top')
#             else:
#                 cage.prepare_initial_topology()
#
#             cage.extract_unique_metal_sites()
#
#             for idx_site in range(len(cage.sites)):
#                 cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#                 cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#                 cage.sites[idx_site].load_fingerprint()
#                 cage.sites[idx_site].set_cutoff()
#             cage.unique_sites = []
#             cage.parametrize(out_coord='crystal_metallicious.pdb', out_topol='crystal_metallicious.top')

# #
# #folder_names = ['Co_cage']

# folder_names = ['megacage']
#
# for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
#
#     for idx, folder_name in enumerate(folder_names):
#
#         print('XXXXXX', folder_name, truncation_scheme_string)
#         basis_folder = '10'
#         # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#         os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#         #name = f"site_{folder_name[:2].title():}0"
#         #metal_name = folder_name[:2].title()
#
#         if folder_name == 'megacage':
#             name = f"site_Pd0"
#             metal_name = 'Pd'
#         else:
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#         library_name = 'merz-opc'
#         if folder_name in uff_folders:
#             library_name = 'uff'
#
#         metal_charges = {metal_name: metal_charges_dic[folder_name]}
#         if folder_name == 'Ru_Pd_cage':
#             metal_charges = {'Ru' : 2, 'Pd' : 2}
#
#         if folder_name=='megacage':
#             cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',
#                                             topol='init_topol/noncovalent_complex.top',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#         else:
#             cage = supramolecular_structure(f'crystal.pdb',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#
#
#         # cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', metal_charges=metal_charges,
#         #                                 topol='init_topol/noncovalent_complex.top', vdw_type=library_name,
#         #                                 search_library=False, truncation_scheme=truncation_scheme_string,
#         #                                 covalent_cutoff=2.5)
#
#         if folder_name!="Co_cage":
#             #cage.prepare_initial_topology()
#             cage.extract_unique_metal_sites()
#             for idx_site in range(len(cage.sites)):
#                 cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#                 cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#                 cage.sites[idx_site].load_fingerprint()
#                 cage.sites[idx_site].set_cutoff()
#             cage.unique_sites = []
#         else:
#             cage.extract_unique_metal_sites()
#             for site in cage.unique_sites:
#                 site.fp_coord_file = f"{site.directory:}/template.pdb"
#                 site.fp_topol_file = f"{site.directory:}/template.top"
#             cage.assign_fingerprints()
#
#         cage.parametrize(out_coord=f'crystal_metallicious_{truncation_scheme_string:}.pdb',
#                          out_topol=f'crystal_metallicious_{truncation_scheme_string:}.top')

#
# folder_names = ['Co_cage']
#
# for idx, folder_name in enumerate(folder_names):
#
#     #print('XXXXXX', folder_name, truncation_scheme_string)
#     basis_folder = '10_uff'
#     # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#     os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#     name = f"site_{folder_name[:2].title():}0"
#     metal_name = folder_name[:2].title()
#
#     library_name = 'merz-opc'
#     if folder_name in uff_folders:
#         library_name = 'uff'
#
#     metal_charges = {metal_name: metal_charges_dic[folder_name]}
#     if folder_name == 'Ru_Pd_cage':
#         metal_charges = {'Ru' : 2, 'Pd' : 2}
#
#     cage = supramolecular_structure('init_topol/noncovalent_complex.pdb', metal_charges=metal_charges,
#                                     topol='init_topol/noncovalent_complex.top', vdw_type=library_name,
#                                     search_library=False, truncation_scheme=None,
#                                     covalent_cutoff=2.5)
#
#     if folder_name!="Co_cage":
#         #cage.prepare_initial_topology()
#         cage.extract_unique_metal_sites()
#         for idx_site in range(len(cage.sites)):
#             cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#             cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#             cage.sites[idx_site].load_fingerprint()
#             cage.sites[idx_site].set_cutoff()
#         cage.unique_sites = []
#     else:
#         cage.extract_unique_metal_sites()
#         for site in cage.unique_sites:
#             site.fp_coord_file = f"{site.directory:}/template.pdb"
#             site.fp_topol_file = f"{site.directory:}/template.top"
#         cage.assign_fingerprints()
#
#     cage.parametrize(out_coord=f'crystal_metallicious.pdb',
#                      out_topol=f'crystal_metallicious.top')
#



#for basis_folder in ['10', 'TZVPP', 'wb97']:
#basis_folder = 'wb97'
#
#for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
#for truncation_scheme_string in ['bond']:
# folder_names= ['Ru_Pd_cage']

# folder_names = ['megacage']
#
# for basis_folder in ['10']:
#     for idx, _ in list(enumerate(folder_names)):
#         folder_name = folder_names[idx]
#         print(folder_name)
#
#         os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}_torsion/{folder_name}/")
#
#         # name = f"site_{folder_name[:2].title():}0"
#         # metal_name = folder_name[:2].title()
#
#         if folder_name == 'megacage':
#             name = f"site_Pd0"
#             metal_name = 'Pd'
#         else:
#             name = f"site_{folder_name[:2].title():}0"
#             metal_name = folder_name[:2].title()
#
#         library_name = 'merz-opc'
#         if folder_name in uff_folders:
#             library_name = 'uff'
#
#         metal_charges = {metal_name: metal_charges_dic[folder_name]}
#         if folder_name == 'Ru_Pd_cage':
#             metal_charges = {'Ru' : 2, 'Pd' : 2}
#
#         # cage = supramolecular_structure(f'crystal.pdb',
#         #                                 metal_charges=metal_charges, vdw_type = library_name, search_library=False)
#         #
#         #cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')
#
#
#         if folder_name=='megacage':
#             cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',
#                                             topol='init_topol/noncovalent_complex.top',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#         else:
#             cage = supramolecular_structure(f'crystal.pdb',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#             cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')
#
#         if folder_name!="Co_cage":
#             cage.extract_unique_metal_sites()
#             for idx_site in range(len(cage.sites)):
#                 cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#                 cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#                 cage.sites[idx_site].load_fingerprint()
#                 cage.sites[idx_site].set_cutoff()
#             cage.unique_sites = []
#         else:
#             cage.extract_unique_metal_sites()
#             for site in cage.unique_sites:
#                 site.fp_coord_file = f"{site.directory:}/template.pdb"
#                 site.fp_topol_file = f"{site.directory:}/template.top"
#             cage.assign_fingerprints()
#
#
#         cage.parametrize(out_coord='crystal_metallicious_impropers.pdb',
#                          out_topol='crystal_metallicious_impropers.top')
#
#
#
#         if folder_name=='megacage':
#             cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',
#                                             topol='init_topol/noncovalent_complex.top',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#         else:
#             cage = supramolecular_structure(f'crystal.pdb',
#                                             metal_charges=metal_charges,
#                                             vdw_type = library_name, search_library=False)
#             cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')
#
#
#
#
#         # cage = supramolecular_structure(f'crystal.pdb', metal_charges=metal_charges,
#         #                                 vdw_type = library_name, search_library=False)
#         #
#         #
#         #
#         #
#         # cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')
#
#         #cage.prepare_initial_topology()
#
#         #cage = supramolecular_structure(f'crystal.pdb',
#         #                                metal_charges=metal_charges,
#         #                                vdw_type = library_name, search_library=False)
#         #cage.prepare_initial_topology()
#
#         if folder_name!="Co_cage":
#             cage.extract_unique_metal_sites()
#             for idx_site in range(len(cage.sites)):
#                 cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#                 cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template_saturated.top'
#                 cage.sites[idx_site].load_fingerprint()
#                 cage.sites[idx_site].set_cutoff()
#             cage.unique_sites = []
#         else:
#             cage.extract_unique_metal_sites()
#             for site in cage.unique_sites:
#                 site.fp_coord_file = f"{site.directory:}/template.pdb"
#                 site.fp_topol_file = f"{site.directory:}/template_saturated.top"
#             cage.assign_fingerprints()
#
#
#         cage.parametrize(out_coord='crystal_metallicious_noimpropers.pdb',
#                          out_topol='crystal_metallicious_noimpropers.top')
#

#
# folder_names = ['Co_mof',  'Zn_mof']
# #folder_names = ['Zn_mof']
#
# metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
#                      'Pd_fujita':2, 'Pd_o':2, 'Ru_Pd_cage':2,   'Zn_knot':2,  'Zn_mof':2,  'Zn_protein':2}
#
# # for basis_folder in ['TZVPP', 'wb97']:
folder_names = ['megacage']

#for basis_folder in ['10_zhang']:#, '10_uff']:

for basis_folder in ['10_uff', '10_zhang']:
    for idx, _ in list(enumerate(folder_names)):
        if not (basis_folder == 'TZVPP' and folder_names[idx] == 'Ru_Pd_cage'):
            folder_name = folder_names[idx]
            print(folder_name)

            os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")

            if folder_name == 'megacage':
                name = f"site_Pd0"
                metal_name = 'Pd'
            else:
                name = f"site_{folder_name[:2].title():}0"
                metal_name = folder_name[:2].title()

            # name = f"site_{folder_name[:2].title():}0"
            # metal_name = folder_name[:2].title()

            library_name = 'merz-opc'
            if folder_name in uff_folders:
                library_name = 'uff'

            metal_charges = {metal_name: metal_charges_dic[folder_name]}

            if folder_name[:2] == 'Zn':
                name_for = {'Zn': (2, 1)}
            elif folder_name[:2] == 'Co':
                name_for = {'Co': (2, 4)}


            #cage = supramolecular_structure('topol/enormous.gro', name_for, topol='topol/enormous.top',
            #                                vdw_type='merz-opc',  search_library=False)

            if folder_name == 'megacage':
                cage = supramolecular_structure('init_topol/noncovalent_complex.pdb',
                                                topol='init_topol/noncovalent_complex.top',
                                                metal_charges=metal_charges,
                                                LJ_type=library_name, search_library=False)
            else:
                cage = supramolecular_structure(f'crystal.pdb',
                                                metal_charges=metal_charges,
                                                LJ_type=library_name, search_library=False)
                cage.prepare_initial_topology(homoleptic_ligand_topol='linker0.top')

            # cage.parametrize(out_topol='out_enormous.top', out_coord='out_enormous.pdb')

            cage.extract_unique_metal_sites()

            for idx_site in range(len(cage.sites)):
                cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
                cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
                cage.sites[idx_site].load_fingerprint()
                cage.sites[idx_site].set_cutoff()
            cage.unique_sites = []
            cage.parametrize(out_coord='crystal_metallicious.pdb', out_topol='crystal_metallicious.top')


# metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
#                      'Pd_fujita':2, 'Pd_o':2, 'Ru_Pd_cage':2,   'Zn_knot':2, 0000 'Zn_mof':2,  'Zn_protein':2}
#
# folder_names = ['Co_mof',  'Zn_mof']
# for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
#
#     for idx, folder_name in enumerate(folder_names):
#
#         print('XXXXXX', folder_name, truncation_scheme_string)
#         basis_folder = '10'
#         # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
#         os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")
#
#         name = f"site_{folder_name[:2].title():}0"
#         metal_name = folder_name[:2].title()
#
#         library_name = 'merz-opc'
#         if folder_name in uff_folders:
#             library_name = 'uff'
#
#         metal_charges = {metal_name: metal_charges_dic[folder_name]}
#         if folder_name[:2] == 'Zn':
#             name_for = {'Zn': (2, 1)}
#         elif folder_name[:2] == 'Co':
#             name_for = {'Co': (2, 4)}
#
#         cage = supramolecular_structure('topol/enormous.gro', name_for, topol='topol/enormous.top',
#                                         vdw_type='merz-opc', search_library=False)
#
#
#         cage.extract_unique_metal_sites()
#         for idx_site in range(len(cage.sites)):
#             cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#             cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#             cage.sites[idx_site].load_fingerprint()
#             cage.sites[idx_site].set_cutoff()
#         cage.unique_sites = []
#
#         cage.parametrize(out_coord=f'crystal_metallicious_{truncation_scheme_string:}.pdb',
#                          out_topol=f'crystal_metallicious_{truncation_scheme_string:}.top')
#

folder_names = ['Zn_mof']

metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
                     'Pd_fujita':2, 'Pd_o':2, 'Ru_Pd_cage':2,   'Zn_knot':2,  'Zn_mof':2,  'Zn_protein':2}

#
# for basis_folder in ['10']:
#     for idx, _ in list(enumerate(folder_names)):
#         folder_name = folder_names[idx]
#         print(folder_name)
#         print('XXXXXX', folder_name, "IMPROP")
#
#         os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}_torsion/{folder_name}/")
#
#         name = f"site_{folder_name[:2].title():}0"
#         metal_name = folder_name[:2].title()
#
#         library_name = 'merz-opc'
#         if folder_name in uff_folders:
#             library_name = 'uff'
#
#         metal_charges = {metal_name: metal_charges_dic[folder_name]}
#
#         if folder_name[:2] == 'Zn':
#             name_for = {'Zn': (2, 1)}
#         elif folder_name[:2] == 'Co':
#             name_for = {'Co': (2, 4)}
#
#         cage = supramolecular_structure('topol/enormous.gro', name_for, topol='topol/enormous.top',
#                                         vdw_type='merz-opc', search_library=False)
#
#         cage.extract_unique_metal_sites()
#         # for site in cage.unique_sites:
#         #     site.fp_coord_file = f"{site.directory:}/template.pdb"
#         #     site.fp_topol_file = f"{site.directory:}/template.top"
#         #
#         # cage.assign_fingerprints()
#
#         for idx_site in range(len(cage.sites)):
#             cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#             cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template.top'
#             cage.sites[idx_site].load_fingerprint()
#             cage.sites[idx_site].set_cutoff()
#         cage.unique_sites = []
#
#         cage.parametrize(out_coord='crystal_metallicious_impropers.pdb',
#                          out_topol='crystal_metallicious_impropers.top')
#
#         cage = supramolecular_structure('topol/enormous.gro', name_for, topol='topol/enormous.top',
#                                         vdw_type='merz-opc', search_library=False)
#
#         cage.extract_unique_metal_sites()
#         for idx_site in range(len(cage.sites)):
#             cage.sites[idx_site].fp_coord_file = f'site_{cage.sites[idx_site].metal_name:}0/template.pdb'
#             cage.sites[idx_site].fp_topol_file = f'site_{cage.sites[idx_site].metal_name:}0/template_saturated.top'
#             cage.sites[idx_site].load_fingerprint()
#             cage.sites[idx_site].set_cutoff()
#         cage.unique_sites = []
#
#         cage.parametrize(out_coord='crystal_metallicious_noimpropers.pdb',
#                          out_topol='crystal_metallicious_noimpropers.top')
#
