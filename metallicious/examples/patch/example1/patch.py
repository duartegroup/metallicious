import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from metallicious.main import patcher
except:
    from cgbind2pmd.main import cgbind2pmd

f='start.gro'
linker_topol='ligand.top'
metal='FE'
metal_charge=2
fingerprint='Fe1'
#fingerprint_style='trunc'

cgbind2gmx = patcher()
cgbind2gmx.name_of_binding_side = fingerprint
#cgbind2gmx.fingerprint_style = fingerprint_style

cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
cgbind2gmx.save(output_coords="cage.gro", output_topol="cage.top")

