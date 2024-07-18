from metallicious.load_fingerprint import guess_fingerprint

def test_guessing():
    fp1 = guess_fingerprint('cage.xyz', 0, metal_name = 'Pd', metal_charge=2, vdw_type='merz-opc')#, fingerprint_guess_list=['Pd2d', 'PdB1', 'PdB2'])
    assert fp1 is not None

    fp2 = guess_fingerprint('cage2.xyz', 184, metal_name = 'Pd', metal_charge=2, vdw_type='merz-opc')#, fingerprint_guess_list=['Pd2d', 'PdB1', 'PdB2'])
    assert fp2 is not None

    fp3 = guess_fingerprint('cage3.xyz', 160, metal_name = 'Pd', metal_charge=2, vdw_type='merz-opc')
    assert fp3 is not None



'''
def load_fp_from_file(filename_fp_coord, filename_fp_topol, fp_style=None, ignore_truncation_warning=True):
    Loads topology and coordinates of the fingerprint. If fingerprint not "full", then it will be trunked to the specified fingerprint style:
    - dihdral (or dih or 3-bond) - truncated to atoms within 3 bond length from metal
    - angle (or ang or 2-bond) - truncated to atoms within 2 bond length from metal
    - bond (or 1-bond) - truncated to atoms within 1 bond length from metal

    :param filename_fp_coord: (string) path to the coordination file of template
    :param filename_fp_topol: (string) path to the topology file of template
    :param fp_style: (string) cutting scheme
    :return:
    
def reduce_site_to_fingerprint(cage_filename, metal_index, syst_fingerprint, cutoff=9, guessing=False, covalent_cutoff=3.0, donors=None):
    Reduces the structure to size of the fingerprint

    :param cage_filename: (string) filename of the input structure
    :param metal_index: (int) metal index
    :param syst_fingerprint: (MDAnalysis.Universe) loaded fingerprint
    :param cutoff: (float) the cut-off from metal center to which atoms are selected for further analysis.
    :param guessing: (bool) if False (default) the function reduce_site... will raise Error if it cannot reduce to fingerprint.
    If True, the function will return False (this is used for guessing of the fingerprint)
    :param covalent_cutoff: (float) cut-off used to detrimne the metal-ligand interactions
    :param donors: (list(string)) element names of the atoms with which metal can interact (e.g. ['N', 'O'])
    :return:
    
def find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_filename, metal_index, metal_name, syst_fingerprint,
                                                              cutoff=9, guessing=False, covalent_cutoff=3.0, donors=None):
    Reduces the metal binding site to the size of the fingerprint (but the ordering of atoms is potentially different).
    Then it finds the best mapping of one structures onto the other

    :param cage_filename: (string) coordination input structure
    :param metal_index: (int) index of metal atom in input structure
    :param metal_name: (string) metal name
    :param syst_fingerprint: (MDAnalysis.Universe) template structure
    :param cutoff: (float) cut-off to reducte the input structure
    :param guessing: (bool) if True program assumes that we try to find the template, and does not give error if not found
    :param covalent_cutoff: (float) cut-off for metal-ligand interaction
    :param donors: (list(str)) atoms with which metal can form bond
    :return:


def search_library_for_fp(metal_name, metal_charge, vdw_type, library_path, fingerprint_guess_list):
    In library path it searches for the possible fingerprints, library should consist of files pairs with coordinates
    and topology, named {metal}_{charge}_{vdw_type}.pdb and {metal}_{charge}_{vdw_type}.top

    :param metal_name: (str) metal name
    :param metal_charge: (int) formal metal charge
    :param vdw_type: (str) name of the dataset used for L-J interactions
    :param library_path: (str) directory to the library of templates
    :param fingerprint_guess_list: (list(str)) predefined list of which templates to check
    :return:


def guess_fingerprint(cage_filename, metal_index, metal_charge, metal_name=None,  fingerprint_guess_list=None,
                      m_m_cutoff=10, vdw_type=None, library_path=f"{os.path.dirname(__file__):s}/library",
                      search_library=True, additional_fp_files=None, fp_style=None, rmsd_cutoff=2, donors=None):

    Finds a template from library of templates

    :param cage_filename: (str) coordination input structure
    :param metal_index: (int) index of metal atom in input structure
    :param metal_name: (string) metal name
    :param metal_charge: (int) metal_charge
    :param fingerprint_guess_list: (list(str)) predefined list of which templates to check
    :param m_m_cutoff: (float) distance between the closest metals
    :param vdw_type: (str) name of the dataset used for L-J interactions
    :param library_path: (str) directory to the library of templates
    :param search_library: (bool) if True it will search library_path for templates
    :param additional_fp_files: (list(str)) additional fingerprints (outside the library directory)
    :param fp_style: (string) cutting scheme
    :param rmsd_cutoff: (float) the guess is accepted if RMSD between metal site and template is below this value
    :param donors: (list(str)) atoms with which metal can form bond
    :return:
'''
