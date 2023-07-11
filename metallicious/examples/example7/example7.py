import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from metallicious.main import patcher
    from metallicious.extract_metal_site import extract_metal_structure
except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure

f='fujita.xyz'
linker_topol='ligand.itp'
metal='Pd'
metal_charge=2
fingerprint='Pd3'
fingerprint_style='full'


#extract_metal_structure(f, metal, "out")

cgbind2gmx = patcher()
cgbind2gmx.name_of_binding_side = fingerprint
cgbind2gmx.fingerprint_style = fingerprint_style
cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))
cgbind2gmx.save(output_coords="out.pdb", output_topol="out.top")



#linker_topol='ligand.top'
#metal='Co'
#metal_charge=2
#fingerprint='Co1'
#fingerprint_style='trunc'

#cgbind2gmx = metallicious()
#cgbind2gmx.name_of_binding_side = fingerprint
#cgbind2gmx.fingerprint_style = fingerprint_style

#cgbind2gmx.from_coords(f, linker_topol, metal, int(metal_charge))

