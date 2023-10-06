from metallicious.parametrize_new_sites import supramolecular_structure
import os

folder_names = ['Pd_Lusby','Co_cage',  'Co_mof',  'Fe_cage',  'Fe_knot',   'Ga_cage',  'Pd_fujita',    'Pd_o', 'Ru_Pd_cage',   'Zn_knot',  'Zn_mof',  'Zn_protein']
uff_folders = ['Ga_cage', 'Ru_Pd_cage' ]

metal_charges_dic = {'Pd_Lusby':2,'Co_cage':2, 'Co_mof':2,  'Fe_cage':2,  'Fe_knot':2,   'Ga_cage':3,
                     'Pd_fujita':2, 'Pd_o':2, 'Ru_Pd_cage':2,   'Zn_knot':2,  'Zn_mof':2,  'Zn_protein':2}

#for basis_folder in ['10', 'TZVPP', 'wb97']:
'''
for basis_folder in ['wb97']:
    #for basis_folder in ['10']:
    for idx, _ in list(enumerate(folder_names)):
        folder_name = folder_names[idx]
        print(folder_name)
        #os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
        os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")

        name = f"site_{folder_name[:2].title():}0"
        metal_name = folder_name[:2].title()

        library_name = 'merz-opc'
        if folder_name in uff_folders:
            library_name = 'uff'

        cage = supramolecular_structure(f'{name:}/saturated_template.xyz',
                                        metal_charges={metal_name: metal_charges_dic[folder_name]},
                                        vdw_type = library_name, search_library=False)

        cage.prepare_initial_topology()
        cage.extract_unique_metal_sites()
        cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
        cage.sites[0].fp_topol_file = f'{name:}/template.top'
        cage.sites[0].load_fingerprint()
        cage.sites[0].set_cutoff()
        cage.unique_sites = []
        cage.parametrize(out_coord='saturated_template.pdb', out_topol='saturated_template.top')

'''
#for basis_folder in ['10', 'TZVPP', 'wb97']:

for truncation_scheme_string in ['dihedral', 'angle', 'bond']:
    #for truncation_scheme_string in [ 'bond']:

    # for basis_folder in ['10']:
    for idx, _ in enumerate(folder_names):

        folder_name = folder_names[idx]

        print('XXXXXX', folder_name, truncation_scheme_string)
        basis_folder = '10'
        # os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
        os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")

        name = f"site_{folder_name[:2].title():}0"
        metal_name = folder_name[:2].title()

        library_name = 'merz-opc'
        if folder_name in uff_folders:
            library_name = 'uff'

        cage = supramolecular_structure(f'{name:}/saturated_template.xyz',
                                        metal_charges={metal_name: metal_charges_dic[folder_name]},
                                        vdw_type=library_name, search_library=False, truncation_scheme=truncation_scheme_string, covalent_cutoff=2.5)

        cage.prepare_initial_topology()
        cage.extract_unique_metal_sites()
        cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
        cage.sites[0].fp_topol_file = f'{name:}/template.top'
        cage.sites[0].load_fingerprint()
        cage.sites[0].set_cutoff()
        cage.unique_sites = []
        cage.parametrize(out_coord=f'saturated_template_{truncation_scheme_string:}.pdb', out_topol=f'saturated_template_{truncation_scheme_string:}.top')


'''
#for basis_folder in ['TZVPP', 'wb97']:
for basis_folder in ['wb97']:
    for idx, _ in list(enumerate(folder_names)):

        folder_name = folder_names[idx]
        print(folder_name)
        #os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_TZVPP/{folder_name}/")
        os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_{basis_folder:}/{folder_name}/")

        name = f"site_{folder_name[:2].title():}0"
        metal_name = folder_name[:2].title()

        library_name = 'merz-opc'
        if folder_name in uff_folders:
            library_name = 'uff'

        cage = supramolecular_structure(f'{name:}/saturated_template.xyz',
                                        metal_charges={metal_name: metal_charges_dic[folder_name]},
                                        vdw_type = library_name, search_library=False)

        cage.prepare_initial_topology()
        cage.extract_unique_metal_sites()
        cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
        cage.sites[0].fp_topol_file = f'{name:}/template_noimpropers.top'
        cage.sites[0].load_fingerprint()
        cage.sites[0].set_cutoff()
        cage.unique_sites = []
        cage.parametrize(out_coord='saturated_template_noimpropers.pdb', out_topol='saturated_template_noimpropers.top')
'''
#cage.parametrize(out_coord='saturated_template_bond.pdb', out_topol='saturated_template_bond.top')

folder_names = ['Co_mof',  'Ga_cage', 'Zn_mof',  'Zn_protein']

folder_names_cages = ['Pd_Lusby', 'Co_cage', 'Fe_cage',  'Fe_knot', 'Ga_cage',  'Pd_fujita', 'Pd_o', 'Ru_Pd_cage', 'Zn_knot']



folder_names_noncages = ['Co_mof', 'Zn_mof', 'Zn_protein']

