from metallicious.extract_metal_site import extract_metal_structure, add_hydrogens

import MDAnalysis
import shutil
import networkx as nx



selection_site = MDAnalysis.Universe("unsaturated_fragment_1.xyz").atoms
with_hydrogens, _ = add_hydrogens(selection_site, "", [24, 25])
G_after = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(with_hydrogens.atoms, with_hydrogens.atoms.positions))
number_bonds_after = len(G_after.edges())
assert number_bonds_after==35

selection_site = MDAnalysis.Universe("unsaturated_fragment_2.pdb").atoms
with_hydrogens, _ = add_hydrogens(selection_site, "Pd", [2, 24, 35, 13])
with_hydrogens.atoms.write("temp.xyz")
G_after = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(with_hydrogens.atoms, with_hydrogens.atoms.positions))
number_bonds_after = len(G_after.edges())
assert number_bonds_after==56

selection_site = MDAnalysis.Universe("unsaturated_fragment_3.pdb").atoms
with_hydrogens, _ = add_hydrogens(selection_site, "FE", [6, 1, 49, 54, 25, 30])
G_after = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(with_hydrogens.atoms, with_hydrogens.atoms.positions, vdwradii={"FE": 0.1}))
number_bonds_after = len(G_after.edges())
assert number_bonds_after==87

selection_site = MDAnalysis.Universe("unsaturated_fragment_4.pdb").atoms
with_hydrogens, _ = add_hydrogens(selection_site, "", [18, 35, 46, 7])
with_hydrogens.atoms.write("temp.xyz")
G_after = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(with_hydrogens.atoms, with_hydrogens.atoms.positions))
number_bonds_after = len(G_after.edges())
assert number_bonds_after==35







extract_metal_structure('tall_cage.pdb', 'tall_cage.top', 'Pd', 'out_tall_cage')
syst = MDAnalysis.Universe("out_tall_cage0/saturated_template.xyz")
assert len(syst.atoms) == 81
shutil.rmtree('out_tall_cage0')

extract_metal_structure('arochain_cage.pdb', 'arochain_cage.top', 'Pd', 'out_arochain_cage')
syst = MDAnalysis.Universe("out_arochain_cage0/saturated_template.xyz")
assert len(syst.atoms) == 129
shutil.rmtree('out_arochain_cage0')

#open topol # TODO

extract_metal_structure('pyryne_cage.pdb', 'pyryne_cage.top', 'Pd', 'out_pyryne_cage')
syst = MDAnalysis.Universe("out_pyryne_cage0/saturated_template.xyz")
assert len(syst.atoms) == 113
shutil.rmtree('out_pyryne_cage0')






# def extract_metal_structure(filename, topol_filename, metal_name, output = None, check_uniquness = True,
#                             all_metal_names= None, covalent_cutoff = 3.0, donors = None, closest_neighbhors = 3):
#     It takes an input supramolecular structure, and it tries extract metal binding site
#
#     :param filename:
#     :param metal_name:
#     :param output:
#     :param check_uniquness: it is used for the strain calculations
#     additional_metals -> metals which are present in topology but removed for parametrization of small molecules
#     :return:



'''

def find_metal_indices(cage, metal_name):
    
    Find indices of metals in the MDAnalysis object
    :param cage: (MDAnalysis.Universe) system in which find MDANalysis object
    :param metal_name: name of the metal
    :return: (list(int), int) list of indices and number of metal (len(list))
    


def find_bound_ligands_nx(cage, metal_index, cutoff=7, cutoff_covalent=3.0, closest_neighbhors=3, neighbhor_cutoff=None,
                          donors=None):
    
    Finds bound ligands to metal, assumes that atoms within cutoff_covalent (default 3) are bound to metal.
    Returns list of sub-graphs of bound ligands; it does not cut ligands (so they can be uneven)

    :param closest_neighbhors:
    :param cage:
    :param donors:
    :param metal_index:
    :param cutoff:
    :param cutoff_covalent:
    :return:
    

def find_closest_and_add_rings(metal_index, cage, bound_ligands, selected_closest_atoms, aromaticity):
    
    This procedure find the atoms which are:
    (a) within dihdedral distance from metal (indicated by metal index)
    (b) if directly connected atom to metal is part of aromatic structure, it will save extract whole aromatic system

    :param metal_index: (int) index of metal
    :param cage: (MDAnalysis.Universe) system which includes the metal
    :param bound_ligands: list(list(int)) list of ligands bound to the metal which consist of lists of ligand's atoms indices
    :param selected_closest_atoms: (list(int)) list of closest atoms to metal (donors)
    :param aromaticity: list(bool) list showing if atoms in (cage) are aromatic
    :return: list(int), list(list(int)), list(list(int)): list of atom indecies of extracted structure, (b) the same
    but indecies separated for diffrent ligands, (c) list of linking (unsaturated) atoms
    

def find_closest_and_add_rings_iterate(metal_indices, cage, all_metal_indecies, covalent_cutoff=3.0, donors=None,
                                       closest_neighbhors=3):
    
    Select closest atoms connected to metal.


    :param metal_indices:       metal_indecies
    :param cage:                MDAnalysis input
    :param G_sub_ligands:
    :return:
    

def check_uniqueness(binding_sites_graphs, site_link_atoms, structure, metal_name, rmsd_cutoff=2):
    # we check if the sites exist in the file
    unique_sites = []
    unique_site_link_atoms = []

    logger.info(f"Uniqueness check: Number of binding sites to check: {len(binding_sites_graphs):}")

    for site, site_link_atom in zip(binding_sites_graphs, site_link_atoms):
        exist = False

        for site_fingerprint in unique_sites:
            logger.info(
                f"\tChecking uniqueness of site and fingerprint with atoms:{len(site):}, {len(site_fingerprint):}")
            # print(f"\tChecking uniqueness of site:{len(site):}, {len(site_fingerprint):}")
            # the structures are sorted because MDAnalysis selection atoms sorts them
            # _, rmsd = map_two_structures(site[0], structure[np.sort(site)], structure[np.sort(site_fingerprint)], metal_name)


            # If number of atoms does not agree then add to unique # TODO number of molecules
            if len(structure[site]) == len(structure[site_fingerprint]):
                _, rmsd = map_two_structures(0, structure[site], structure[site_fingerprint],
                                             metal_name)
            else:
                #print("NOT The same")
                rmsd = 1e10

            logger.info(f"RMSD between current structure and already check structures: {rmsd:}")
            # print(f"RMSD between current structure and already check structures: {rmsd:}")
            if rmsd < rmsd_cutoff:
                exist = True
                break

        if exist == False:
            logger.info(f"Adding unique site # {len(unique_sites):}, {[len(a) for a in unique_sites]:}")
            unique_sites.append(site)
            unique_site_link_atoms.append(site_link_atom)

    logger.info(f"Number of unique sites: {len(unique_sites):}")
    return unique_sites, unique_site_link_atoms



def rotation_matrix_from_vectors(vec1, vec2):
    Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix



def add_hydrogens(selection_site, metal_name, add_atoms_to_this_atom):
    
    Adds hydrogen to unsaturated carbon

    :param selection_site:
    :param metal_name:
    :param add_atoms_to_this_atom:
    :return:
    



def find_ligand_pattern(new_syst, ligands_nodes):
    
    Check if the site is composed of the same ligand fragments (we want to symmetrize paramters for the same sites)

    :param new_syst:
    :param ligands_nodes:
    :return:
    


def renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern, extra_atoms,
                    link_atoms, covalent_cutoff=3):
    
    Renumbers atoms in ligands, that the identical one have the same numbering, and they are grouped

    :param new_syst:
    :param metal_name:
    :param ligands_atoms_membership:
    :param unique_ligands:
    :param unique_ligands_pattern:
    :param extra_atoms:
    :param link_atoms:
    :return:
    


def read_and_reoder_topol_and_coord(filename, topol_filename, metal_name, all_metal_names=None):
    if topol_filename is None:
        raise ValueError("[!] Topology file missing!")


def extract_metal_structure(filename, topol_filename, metal_name, output = None, check_uniquness = True,
                            all_metal_names= None, covalent_cutoff = 3.0, donors = None, closest_neighbhors = 3):
    It takes an input supramolecular structure, and it tries extract metal binding site

    :param filename:
    :param metal_name:
    :param output:
    :param check_uniquness: it is used for the strain calculations
    additional_metals -> metals which are present in topology but removed for parametrization of small molecules
    :return:


'''