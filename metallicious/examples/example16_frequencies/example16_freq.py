from metallicious.parametrize_new_sites import supramolecular_structure
import os



folder_names = ['Pd_Lusby','Co_cage',  'Co_Lusby',  'Co_mof',  'Fe_cage',  'Fe_knot',   'Ga_cage',  'Pd_fujita',    'Pd_o', 'Ru_Pd_cage',   'Zn_knot',  'Zn_mof',  'Zn_protein']

uff_folders = ['Ga_cage', 'Ru_Pd_cage' ]


for idx, _  in enumerate(folder_names):
    #idx=0
    folder_name = folder_names[idx]
    os.chdir(f"/u/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/{folder_name}/")

    name = f"site_{folder_name[:2].title():}0"
    metal_name = folder_name[:2].title()

    library_name = 'merz-opc'
    if folder_name in uff_folders:
        library_name = 'uff'

    #cage = supramolecular_structure(f'{name:}/saturated_template.xyz',  metal_charges={metal_name: 2}, vdw_type = library_name,
    #                                    search_library=False, truncation_scheme='dihedral')
    cage = supramolecular_structure(f'{name:}/saturated_template.xyz',  metal_charges={metal_name: 2}, vdw_type = library_name,
                                    search_library=False, truncation_scheme='dihedral')


    cage.prepare_initial_topology()
    cage.extract_unique_metal_sites()
    cage.sites[0].fp_coord_file = f'{name:}/template.pdb'
    cage.sites[0].fp_topol_file = f'{name:}/template.top'
    cage.sites[0].load_fingerprint()
    cage.sites[0].set_cutoff()
    cage.unique_sites = []
    #cage.parametrize(out_coord='saturated_template_dihedral.pdb', out_topol='saturated_template_dihedral.top')
    cage.parametrize(out_coord='saturated_template.pdb', out_topol='saturated_template.top')




# os.chdir('/home/fd05/fd/chem1540/Research/2021_11_02_Pullen/classic_cages_07_10/Zn_mof')
# cage = supramolecular_structure('topol/enormous.gro', {'Zn': (2,1)}, topol='topol/enormous.top', vdw_type='merz-opc')
# cage.parametrize(out_topol='out_enormous.top', out_coord='out_enormous.pdb')
