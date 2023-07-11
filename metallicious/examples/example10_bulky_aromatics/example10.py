import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from metallicious.main import patcher
    from metallicious.extract_metal_site import extract_metal_structure
except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure

f='cage3.xyz'

metal='Pd'
metal_charge=2
fingerprint_style='trunc'
extract_metal_structure(f, metal, "out")

