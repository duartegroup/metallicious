import parmed as pmd
import MDAnalysis


def compare_topology_and_coords(topol_filename, coord_filename):
    syst = MDAnalysis.Universe(coord_filename)
    topol = pmd.load_file(topol_filename)

    assert len(syst.atoms) == len(topol.atoms)



def check_if_orca_available():
    try:
        import autode as ade
    except:
        raise ImportError("No autode found")

    method = ade.methods.ORCA()

    if method.is_available is False:
        raise NameError("For parametrization of templates, QM software ORCA is required")


def check_if_psriresp_available():
    try:
        import psiresp
    except:
        raise ImportError("No psiresp found")

def check_if_parametrization_modules_available():
    check_if_orca_available()
    check_if_psriresp_available()