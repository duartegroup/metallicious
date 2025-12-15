'''
Parametrization of the cage, using non-bonded topology of the whole cage

Equivalent bash command:
metallicious -f ru_pd.pdb -p ru_pd.top -vdw_type uff -metal_and_charges Pd 2 Ru 2
'''

from metallicious import supramolecular_structure

cage = supramolecular_structure(f'start.pdb',
                                topol = 'start.top',
                                metal_charges={'Cu': 2},
                                LJ_type = 'merz-opc',
                                search_library=False,
                                donors=['N', 'O', 'Cu'],)

cage.extract_unique_metal_sites() # extracts template structures (needed to detect sites)
cage.sites[0].fp_coord_file = 'template.pdb' # loads coordinates of the template
cage.sites[0].fp_topol_file = 'template.top' # loads force-field parameters of the template
cage.sites[0].load_fingerprint() # loads template into the class
cage.sites[0].set_cutoff() # reads the radius of the template
cage.unique_sites = [] # overwrites list with templates to parameterize (no template parametrization needed)

cage.parametrize(out_coord='small_site.pdb', out_topol='small_site.top')
