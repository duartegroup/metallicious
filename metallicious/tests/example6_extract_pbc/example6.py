'''import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/metallicious/')


try:
    from main import metallicious
    from extract_metal_site import extract_metal_structure
except:
    from metallicious.main import metallicious
    from metallicious.extract_metal_site import extract_metal_structure

f='Zn_tet.pdb'
metal='Zn'
extract_metal_structure(f, metal, "out")

'''

from metallicious.parametrize_new_sites import supramolecular_structure


import os

cage = supramolecular_structure('noncovalent_complex.pdb', {'Co': (2, 1)}, topol='noncovalent_complex.pdb', vdw_type='merz-opc')

cage.extract_unique_metal_sites()

site = cage.unique_sites[0]
os.chdir(site.directory)
print("---------------------------------------------------------------------------")
site.create_initial_topol()



#site = cage.sites[0]
#os.chdir(site.filename)
##site.create_initial_topol()
#site.seminario()

#site.partial_charge()
#site.cut_extra_atoms_and_save()

print("")
