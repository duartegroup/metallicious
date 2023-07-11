try:
    from metallicious.main import patcher
except:
    from cgbind2pmd.main import cgbind2pmd

f='square.xyz'
linker_topol='linker.top'
metal='Pd'
metal_charge=2
fingerprint='Pd3'
fingerprint_style='angle'


cgbind2gmx = patcher()
cgbind2gmx.name_of_binding_side = fingerprint
cgbind2gmx.fingerprint_style = fingerprint_style
cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
cgbind2gmx.save(output_coords="out.pdb", output_topol="out.top")


