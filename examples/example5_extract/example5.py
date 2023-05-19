import sys
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

