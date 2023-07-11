import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from metallicious.main import patcher
except:
    from cgbind2pmd.main import cgbind2pmd

f='Ag_tet.pdb'
linker_topol='linker_topol.top'
metal='Ag'
metal_charge=1
fingerprint='Ag1'
fingerprint_style='trunc'

cgbind2gmx = patcher()
cgbind2gmx.name_of_binding_side = fingerprint
cgbind2gmx.fingerprint_style = fingerprint_style

cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
cgbind2gmx.save(output_coords="out.pdb", output_topol="out.top")
