import os
import parmed as pmd
import MDAnalysis
import numpy as np
import networkx as nx

from MDAnalysis.lib.distances import distance_array

from metallicious.extract_metal_site import find_bound_ligands_nx
from metallicious.log import logger
from metallicious.mapping import map_two_structures
from metallicious.utils import strip_numbers_from_atom_names

def load_fp_from_file(filename_fp_coord, filename_fp_topol, fp_style=None, ignore_truncation_warning=True):
    '''
    Loads topology and coordinates of the fingerprint. If fingerprint not "full", then it will be trunked to the specified fingerprint style:
    - dihdral (or dih or 3-bond) - truncated to atoms within 3 bond length from metal
    - angle (or ang or 2-bond) - truncated to atoms within 2 bond length from metal
    - bond (or 1-bond) - truncated to atoms within 1 bond length from metal

    :param filename_fp_coord: (string) path to the coordination file of template
    :param filename_fp_topol: (string) path to the topology file of template
    :param fp_style: (string) cutting scheme
    :return:
    '''

    topol = pmd.load_file(filename_fp_topol)
    syst_fingerprint = MDAnalysis.Universe(filename_fp_coord)

    residual_charges = np.array([np.sum([atom.charge for atom in residue]) for residue in topol.residues])

    if fp_style is None:
        # do not truncate
        return topol, syst_fingerprint

    elif fp_style in ['dihedral', 'dih', 'angle', 'ang', 'bond']:
        bonds = [(bond.atom1.idx, bond.atom2.idx) for bond in topol.bonds]
        G = nx.Graph(bonds)

        if fp_style == 'dihedral' or fp_style == 'dih' or fp_style == '3-bond':
            atoms_bound_to_metal_by_bonded = list(nx.generators.ego_graph(G, 0, radius=3).nodes)

        elif fp_style == 'angle' or fp_style == 'ang' or fp_style == '2-bond':
            atoms_bound_to_metal_by_bonded = list(nx.generators.ego_graph(G, 0, radius=2).nodes)

        elif fp_style == 'bond' or fp_style == '1-bond':
            atoms_bound_to_metal_by_bonded = list(nx.generators.ego_graph(G, 0, radius=1).nodes)
            # atoms_bound_to_metal_by_bonded = list(
            #     set(np.concatenate([[bond.atom1.idx, bond.atom2.idx] for bond in metal_topology.bonds])))

        atoms_to_strip = []
        residual_charge = 0.0
        for atom in topol.atoms:
            if atom.idx not in atoms_bound_to_metal_by_bonded:
                atoms_to_strip.append(atom.idx)
                residual_charge += atom.charge

        residual_charge /= len(atoms_bound_to_metal_by_bonded)

        for idx in atoms_bound_to_metal_by_bonded:
            topol.atoms[idx].charge += residual_charge

        topol.strip(f"@{','.join(list(map(str, np.array(atoms_to_strip) + 1))):s}")

        residual_charges_cut = np.array([np.sum([atom.charge for atom in residue]) for residue in topol.residues])
        std = np.std(residual_charges_cut - residual_charges)
        logger.info(f"The truncation resulted in error of partial charges: {std:.02f}")

        if std > 0.3 and ignore_truncation_warning==False:
            raise ValueError(f"WARNING the truncation results in large overall change of charge of residues: changing truncation scheme might help; std: {std:}")


        new_syst_fingerprint = syst_fingerprint.select_atoms(
            f"not index {' '.join(list(map(str, np.array(atoms_to_strip)))):s}")
        new_syst_fingerprint.write("temp.pdb")
        new_syst_fingerprint = MDAnalysis.Universe("temp.pdb")



        return topol, new_syst_fingerprint.atoms
    else:
        raise ValueError("Incorrect value of truncation scheme")


def reduce_site_to_fingerprint(cage_filename, metal_index, syst_fingerprint, cutoff=9, guessing=False, covalent_cutoff=3.0, donors=None):
    '''
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
    '''


    cage = MDAnalysis.Universe(cage_filename)

    metal_type = syst_fingerprint.atoms[0].type
    bonds_with_metal = MDAnalysis.topology.guessers.guess_bonds(syst_fingerprint.atoms,
                                                                syst_fingerprint.atoms.positions,
                                                                vdwradii={metal_type: covalent_cutoff})

    G_fingerprint_with_metal = nx.Graph(bonds_with_metal)
    neighbhor_cutoff = nx.eccentricity(G_fingerprint_with_metal, v=0)

    syst_fingerprint_no_metal = syst_fingerprint.atoms[1:]

    bonds = MDAnalysis.topology.guessers.guess_bonds(syst_fingerprint_no_metal.atoms,
                                                     syst_fingerprint_no_metal.atoms.positions)
    G_fingerprint = nx.Graph()

    if len(bonds) > 0:  # fingerprint with ligands larger than one atom
        G_fingerprint = nx.Graph(bonds)
    elif len(bonds) == 0:  # not bonds, this has to be minimal fingerprint, only donor atoms
        G_fingerprint.add_nodes_from(syst_fingerprint_no_metal.atoms.indices)

    nx.set_node_attributes(G_fingerprint, {atom.index: atom.name[0] for atom in syst_fingerprint_no_metal.atoms},
                           "name")
    G_fingerprint_subs = [G_fingerprint.subgraph(a) for a in nx.connected_components(G_fingerprint)]
    logger.info(f"\t\t\t[ ] Mapping fingerprint to metal center: {metal_index:d}")

    G_sub_cages, closest_atoms = find_bound_ligands_nx(cage, metal_index, cutoff=cutoff, neighbhor_cutoff=neighbhor_cutoff, cutoff_covalent=covalent_cutoff, donors=donors)


    closest_atoms = np.concatenate(closest_atoms)
    number_ligands_bound = len(G_sub_cages)

    if len(G_sub_cages) != len(G_fingerprint_subs) and guessing:
        logger.info(f"\t\t\t[!] Not the same number of sites {guessing:}, structure: {number_ligands_bound:d} vs. fingerprint {len(G_fingerprint_subs):d}")
        return False

    elif len(G_sub_cages) != len(G_fingerprint_subs) and not guessing:
        if len(G_sub_cages)>len(G_fingerprint_subs):
            raise ValueError(
                f"[!] Not the same number of sites {len(G_sub_cages):}!={len(G_fingerprint_subs):} guessing={guessing:}. Probably, some of the hydrogen/carbons are too close. Try increasing cutoff")

        raise ValueError(f"[!] Not the same number of sites {len(G_sub_cages):}!={len(G_fingerprint_subs):} guessing={guessing:}")

    selected_atoms = []
    end_atoms = []

    fingerprint_ids = []

    for G_sub_cage in G_sub_cages:
        largest_common_subgraph = []
        finerprint_idx = None
        trial = 0
        for G_idx, G_fingerprint_sub in enumerate(G_fingerprint_subs):

            # ISMAGS is a slow algorithm, so we want to makes sure that we match correct substructure, we make simple assesments:
            # Let's try to use rdkit:


            mol1 = MDAnalysis.Merge(cage.select_atoms(f'index {" ".join(list(map(str, list(G_sub_cage)))):s}'))
            if cage.dimensions is not None: # this is only for cage neccessary, as templates should not be crossing pbc
                mol1.dimensions = cage.dimensions


            #mol1 = cage.select_atoms(f'index {" ".join(list(map(str, list(G_sub_cage)))):s}')
            if not hasattr(mol1.atoms[0], 'element'):
                guessed_elements = MDAnalysis.topology.guessers.guess_types(strip_numbers_from_atom_names(mol1.atoms.names))
                mol1.add_TopologyAttr('elements', guessed_elements)
            # Rremove dummy atoms (i.g., OPC, TIP4P) which break rdkit
            mol1 = mol1.select_atoms("not element DUMMY")


            mol2 = MDAnalysis.Merge(syst_fingerprint_no_metal.select_atoms(f'index {" ".join(list(map(str, list(G_fingerprint_sub)))):s}'))

            cage_sub_rdkit = mol1.atoms.convert_to("RDKIT", NoImplicit=False) # NoImplicit is use to allow implicit hydrogens (but if they are there it does not remove them!). Important side effect is that it does not bond orders etc. allowing for comparison of structures

            fp_sub_rdkit = mol2.atoms.convert_to("RDKIT", NoImplicit=False)
            passed_conditions = cage_sub_rdkit.HasSubstructMatch(fp_sub_rdkit)

            if passed_conditions:
                ismags = nx.isomorphism.ISMAGS(G_sub_cage, G_fingerprint_sub,
                                               node_match=lambda n1, n2: n1['name'] == n2['name'])

                for largest_common_subgraph_iter in ismags.largest_common_subgraph(symmetry=False):
                    if len(largest_common_subgraph_iter) < len(G_fingerprint_sub):
                        break

                    is_donor_included = [closest_atom in largest_common_subgraph_iter for closest_atom in
                                         closest_atoms]

                    if any(is_donor_included):
                        if len(largest_common_subgraph_iter) > len(largest_common_subgraph):
                            largest_common_subgraph = largest_common_subgraph_iter
                            finerprint_idx = G_idx
                            logger.debug(f"Found pattern which has all donor atoms {largest_common_subgraph:}")
                            break
        # else:
        #    print("Mismatch by size")
        # else:
        if finerprint_idx is None:
            trial += 1
            # we make sure that we find the matching ligand

            if trial == len(G_fingerprint_sub):
                if guessing:
                    return False
                else:
                    assert trial < len(G_fingerprint_sub)

        else:
            fingerprint_ids.append(finerprint_idx)

        if guessing:  # if we want to guess site, that might be not fulfilled, and that means it is not the correct site
            if len(largest_common_subgraph) == 0:
                return False
            elif not len(G_fingerprint_subs[finerprint_idx]) == len(largest_common_subgraph):
                return False
        else:
            # Check if the subgraph is the same size as fingerprint,
            # print(largest_common_subgraph)

            assert len(largest_common_subgraph) > 0
            assert len(G_fingerprint_subs[finerprint_idx]) == len(largest_common_subgraph)

        selected_atoms += largest_common_subgraph.keys()

        # Let's find the end atoms, they have different degree: in cage they are connected, in the template not
        G_sub_cage_degree = G_sub_cage.degree
        G_fingerprint_subs_degree = G_fingerprint_subs[finerprint_idx].degree

        # We cut the last atom
        for key in largest_common_subgraph:
            if G_sub_cage_degree[key] != G_fingerprint_subs_degree[largest_common_subgraph[key]]:
                end_atoms.append(key) #

    connected_cut_system = cage.select_atoms(f"index {metal_index:}") + cage.select_atoms(
        f"index {' '.join(map(str, selected_atoms)):}")

    return connected_cut_system #, end_atoms


def find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_filename, metal_index, metal_name, syst_fingerprint,
                                                              cutoff=9, guessing=False, covalent_cutoff=3.0, donors=None):
    '''
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
    '''

    connected_cut_system = reduce_site_to_fingerprint(cage_filename, metal_index, syst_fingerprint, cutoff=cutoff, guessing=guessing, covalent_cutoff=covalent_cutoff, donors=donors)

    if connected_cut_system:  # the results are legit, we take them as input
        best_mapping, best_rmsd = map_two_structures(metal_index, connected_cut_system, syst_fingerprint, metal_name=metal_name)

    else:  # we are guessing, and we missed
        return None, 1e10

    return best_mapping, best_rmsd


def search_library_for_fp(metal_name, metal_charge, vdw_type, library_path, fingerprint_guess_list):
    '''
    In library path it searches for the possible fingerprints, library should consist of files pairs with coordinates
    and topology, named {metal}_{charge}_{vdw_type}.pdb and {metal}_{charge}_{vdw_type}.top

    :param metal_name: (str) metal name
    :param metal_charge: (int) formal metal charge
    :param vdw_type: (str) name of the dataset used for L-J interactions
    :param library_path: (str) directory to the library of templates
    :param fingerprint_guess_list: (list(str)) predefined list of which templates to check
    :return:
    '''
    # Check available sites in the library directory
    all_fingerprints_names = []
    for file in os.listdir(library_path):
        if file.endswith('.pdb'):
            all_fingerprints_names.append(file[:-4])

    # We choose only the one which have the same name
    fingerprints_names = {}
    if fingerprint_guess_list is not None:
        for name in fingerprint_guess_list:
            if name not in all_fingerprints_names:
                print("Fingerprint not available")
                raise ValueError(f"Selected fingerprint ({name:}) not avaialable")
            else:
                fingerprints_names.append(name)
    else:
        for name in all_fingerprints_names:
            if len(name.split('_')) > 2:
                if metal_name.title() in name:
                    if metal_charge == int(name.split('_')[1]):
                        if vdw_type is None:
                            fingerprints_names[
                                name] = f"{library_path:}/{name:}.pdb"
                        elif vdw_type==name.split('_')[2]:  # if vdw_type is specified we only search for specific files
                            fingerprints_names[name] = f"{library_path:}/{name:}.pdb"
            else:
                logger.warning(f"Fingerprint file name ({name}) in incorrect format")

    return fingerprints_names


def guess_fingerprint(cage_filename, metal_index, metal_charge, metal_name=None,  fingerprint_guess_list=None,
                      m_m_cutoff=10, vdw_type=None, library_path=f"{os.path.dirname(__file__):s}/library",
                      search_library=True, additional_fp_files=None, fp_style=None, rmsd_cutoff=2, donors=None):
    '''
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

    logger.info(f"Guessing template for {metal_name}{metal_charge}+[{metal_index}]")

    fp_files = {}
    if search_library:
        fp_files = {**fp_files,
                    **search_library_for_fp(metal_name, metal_charge, vdw_type, library_path, fingerprint_guess_list)}

    if additional_fp_files is not None:
        fp_files = {**additional_fp_files, **fp_files}

    if len(fp_files) == 0:
        logger.info(f"Not found templates with name including {metal_name:} with vdw type: {vdw_type:}")
        return False

    logger.info(f"\tSelected template files: {fp_files.keys():}")
    rmsd_best = 1e10
    name_of_binding_side = None
    for finerprint_name in fp_files:
        logger.info(f"\t\t[ ] Guessing fingerprint {finerprint_name:s}")

        #syst_fingerprint = MDAnalysis.Universe()
        _, syst_fingerprint = load_fp_from_file(f'{fp_files[finerprint_name]}', f'{fp_files[finerprint_name].replace(".pdb", ".top")}', fp_style=fp_style)

        # we need to find what is smaller
        metal_position = syst_fingerprint.atoms[0].position
        nometal_position = syst_fingerprint.atoms[1:].positions
        cutoff = np.min([np.max(distance_array(metal_position, nometal_position)) + 2.0, m_m_cutoff])
        # cutoff = 0.5*(cutoff+self.m_m_cutoff) # we make it 75% close to metal

        _, rmsd = find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_filename, metal_index, metal_name,
                                                                            syst_fingerprint, guessing=True,
                                                                            cutoff=cutoff, donors=donors)

        if rmsd < rmsd_best:
            rmsd_best = rmsd
            name_of_binding_side = fp_files[finerprint_name].replace(".pdb", "")
            # self.ligand_cutoff = cutoff
        logger.info(f"\t\t\t[ ] RMSD {rmsd:f}")

    if name_of_binding_side is None:
        logger.info(f"\t[-] Not found any fingerprints for {metal_name:} with vdw type: {vdw_type:}")
        return False

    logger.info(f"\t[+] Best fingerprint {name_of_binding_side:s} rmsd: {rmsd_best:f}")
    if rmsd_best > rmsd_cutoff:
        logger.info(f"\t[!] Rmsd of best fingerprint is above {rmsd_cutoff:}, figerprint not guessed")
        return False
    return name_of_binding_side