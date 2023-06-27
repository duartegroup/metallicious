'''import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from main import cgbind2pmd
    from extract_metal_site import extract_metal_structure
except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure

f='Zn_tet.pdb'
metal='Zn'
extract_metal_structure(f, metal, "out")

'''
import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


from parametrize_new_sites import supramolecular_structure


import os

cage = supramolecular_structure('Co_tet.pdb', {'Co': (2,1)}, vdw_type='merz-tip3p')
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
