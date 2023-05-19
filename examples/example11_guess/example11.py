
try:
    from main import cgbind2pmd
except:
    from cgbind2pmd.main import cgbind2pmd

#f='start.pdb'

f='cage.xyz'
#f='template.pdb'
#f='crystal.pdb'
#f='cage_m6l12_gaff.pdb'
#linker_topol='new_topol.top'
linker_topol=None
#linker_topol='new_topol.top'

metal='Pd'
metal_charge=2
fingerprint='Pd2d'
fingerprint_style='dih'

cgbind2gmx = cgbind2pmd()
cgbind2gmx.name_of_binding_side = fingerprint
cgbind2gmx.fingerprint_style = fingerprint_style

cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
cgbind2gmx.save(output_coords="out.pdb", output_topol="out.top")


## TODO:
# cgbind2pmd -f old_cage.xyz -linker_topol ../temp/linker.itp -metal PD -metal_charge 2 -fingerprint Pd2
