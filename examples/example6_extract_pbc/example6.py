import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from main import cgbind2pmd
    from extract_metal_site import extract_metal_structure
except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure

f='Zn_tet.pdb'
linker_topol='linker.top'
metal='Zn'
metal_charge=2
#fingerprint='Zn'
#fingerprint_style='trunc'


extract_metal_structure(f, metal, "out")


#cgbind2gmx.name_of_binding_side = fingerprint
#cgbind2gmx.fingerprint_style = fingerprint_style

#cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
