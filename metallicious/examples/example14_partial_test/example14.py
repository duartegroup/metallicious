import sys
sys.path.insert(0, '/home/fd05/fd/chem1540/github/cgbind2pmd/')


try:
    from metallicious.main import patcher
    from metallicious.extract_metal_site import extract_metal_structure
    from metallicious.initial_site import create_initial_topol
except:
    from cgbind2pmd.main import cgbind2pmd
    from cgbind2pmd.extract_metal_site import extract_metal_structure
    from cgbind2pmd.initial_site import create_initial_topol


#metal_name = 'Pd'
#metal_charge = 2
#vdw_type ='merz-tip3p'

#topols = create_initial_topol(metal_name, metal_charge, vdw_data_name=vdw_type)