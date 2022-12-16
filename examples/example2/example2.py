
try:
    from main import cgbind2pmd
except:
    from cgbind2pmd.main import cgbind2pmd

f='new_cage.gro'
linker_topol='linker.itp'
metal='Pd'
metal_charge=2
fingerprint='Pd2'
#fingerprint_style='trunc'

cgbind2gmx = cgbind2pmd()
cgbind2gmx.name_of_binding_side = fingerprint
#cgbind2gmx.fingerprint_style = fingerprint_style

cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))


# TODO:
# cgbind2pmd -f old_cage.xyz -linker_topol ../temp/linker.itp -metal PD -metal_charge 2 -fingerprint Pd2