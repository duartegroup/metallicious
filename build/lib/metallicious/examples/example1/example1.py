try:
    from metallicious.main import cgbind2pmd
except:
    from cgbind2pmd.main import cgbind2pmd

f='cage_m6l12_ob.xyz'
linker_topol='linker.top'
metal='Pd'
metal_charge=2
fingerprint='Pd1'
fingerprint_style='trunc'


cgbind2gmx = cgbind2pmd()
cgbind2gmx.name_of_binding_side = fingerprint
cgbind2gmx.fingerprint_style = fingerprint_style
cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
