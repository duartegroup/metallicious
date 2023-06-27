import os

import networkx as nx
import MDAnalysis
from itertools import permutations
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from networkx.algorithms import isomorphism
import numpy as np
import re

import MDAnalysis.transformations


try:
    from cgbind2pmd.log import logger
except:
    from log import logger


def syst_to_graph(atoms, vdwradii):
    bonds = MDAnalysis.topology.guessers.guess_bonds(atoms, atoms.positions)
    G_fingerprint = nx.Graph()
    if len(bonds) > 0:  # fingerprint with ligands larger than one atom
        G_fingerprint = nx.Graph(bonds, vdwradii=vdwradii)
    elif len(bonds) == 0:  # not bonds, this has to be minimal fingerprint, only donor atoms
        G_fingerprint.add_nodes_from(atoms.indices)

    nx.set_node_attributes(G_fingerprint, {atom.index: atom.name[0] for atom in atoms.atoms}, "name")
    return G_fingerprint


def unwrap(syst, metal_type):
    '''
    This tries to unwrap the trajectory through PBC
    
    :param syst:
    :return:
    '''
    new_universe = syst.universe.copy()

    G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(syst.atoms,
                                                               syst.atoms.positions,
                                                               vdwradii={metal_type: 3},
                                                               box=new_universe.dimensions))
    bonds = list(G_cage.edges)
    new_universe.add_TopologyAttr('bonds', bonds)
    new_syst = new_universe.atoms[syst.atoms.indices]

    transform = MDAnalysis.transformations.unwrap(new_syst)
    new_universe.trajectory.add_transformations(transform)

    return new_syst

def int_list_to_str(lista):
    return " ".join(list(map(str, lista)))
    

def map_two_structures(metal_index, connected_cut_system, syst_fingerprint, metal_name, end_atoms=[], periodic=False): #TODO why need metal_index and metal_name (?)
    # usually fingerprint suppose to not have the PBC, but for quesing purposes it can have

    syst_fingerprint_pbc = syst_fingerprint
    if syst_fingerprint.universe.dimensions is not None:
        syst_fingerprint_pbc = unwrap(syst_fingerprint, syst_fingerprint[0].type)

    # select_atoms automatically tries to sort things, but we assume that metal is first
    syst_fingerprint_heavy_atoms = syst_fingerprint_pbc.select_atoms(f"not name H*", sorted=False)

    G_fingerprint_heavy_atoms = syst_to_graph(syst_fingerprint_heavy_atoms.atoms[1:], vdwradii={metal_name: 1.0})

    G_fingerprint_subs_heavy_atoms = [G_fingerprint_heavy_atoms.subgraph(a) for a in
                                      nx.connected_components(G_fingerprint_heavy_atoms)]

    connected_cut_system_pbc = connected_cut_system
    if connected_cut_system.universe.dimensions is not None:
        connected_cut_system_pbc = unwrap(connected_cut_system, connected_cut_system[0].type)

    connected_cut_system_pbc_heavy_atoms = connected_cut_system_pbc.select_atoms("not name H*", sorted=False)

    # no_metal_heavy_atoms = connected_cut_system_pbc_heavy_atoms.select_atoms(f"not index {metal_index:d}") # TODO remove in future
    no_metal_heavy_atoms = connected_cut_system_pbc_heavy_atoms[1:]
    
    G_site_heavy_atoms = syst_to_graph(no_metal_heavy_atoms, vdwradii={metal_name: 1.0})

    # G_site_heavy_atoms = nx.Graph(
    #    MDAnalysis.topology.guessers.guess_bonds(no_metal_heavy_atoms.atoms, no_metal_heavy_atoms.atoms.positions))
    # nx.set_node_attributes(G_site_heavy_atoms, {atom.index: atom.name[0] for atom in no_metal_heavy_atoms.atoms}, "name")
    G_site_subs_heavy_atoms = [G_site_heavy_atoms.subgraph(a) for a in nx.connected_components(G_site_heavy_atoms)]

    number_ligands_bound = len(G_site_subs_heavy_atoms)

    best_rmsd = 1e10
    best_mapping = None

    indecies2number = {}
    for a, b in enumerate(connected_cut_system_pbc_heavy_atoms.indices):
        indecies2number[b] = a

    # We find which atom is the metal in new site TODO remove, this is not necessary, we assume is the first one
    # new_metal_index = np.where(connected_cut_system_pbc_heavy_atoms.indices == metal_index)[0][0]
    permutations_list = permutations(range(number_ligands_bound))


    # We check all permutation of residue number
    for perm in list(permutations_list):
        mappings = []
        # some metal sites are not symmetric, if they are not we don't check them
        length_match = True

        # we check all permutation of atom numbers in residues
        for lig_a, lig_b in enumerate(perm):
            # if length of ligands is different, then this is definitely wrong permutation
            if len(G_fingerprint_subs_heavy_atoms[lig_a]) != len(G_site_subs_heavy_atoms[lig_b]):
                mappings = []
                length_match = False
                break

            iso = isomorphism.GraphMatcher(G_fingerprint_subs_heavy_atoms[lig_a], G_site_subs_heavy_atoms[lig_b],
                                           node_match=lambda n1, n2: n1['name'] == n2['name'])
            mappings.append([subgraph_mapping for subgraph_mapping in iso.subgraph_isomorphisms_iter()])

        if length_match:
            # we tabulate all the permutations (of residues) of permutations (of atoms)
            all_possible_mappings = []

            def recursive(mapping, index, new_mapping):
                if index == len(mapping):
                    all_possible_mappings.append(new_mapping)
                else:
                    for map_1 in mapping[index]:
                        copy = new_mapping.copy()
                        copy.append(map_1)
                        recursive(mapping, index + 1, copy)

            recursive(mappings, 0, [])
            #logger.info(f"Permutations to check:{len(all_possible_mappings):}")
        
            for mapping_idx in range(len(all_possible_mappings)): # TODO change range(len to enumerate
                # concatenate the mapping into one dictionary
                concatenated_mapping = {}
                for ligand in all_possible_mappings[mapping_idx]:
                    concatenated_mapping.update(ligand)

                # Reorder left hand side
                reversed_concatenated_mapping = {}
                for d in sorted(concatenated_mapping.keys()):
                    reversed_concatenated_mapping[d] = concatenated_mapping[d]

                reordered_system = connected_cut_system_pbc_heavy_atoms.atoms[
                    [0] + [indecies2number[value] for value in reversed_concatenated_mapping.values()]]  # Find where metal is
                # make sure that the type is preseved (esptially challening with names starting with C and N)
                reordered_system.atoms[0].type = metal_name
                reordered_system.atoms[0].mass = 0.0

                # MDAnalysis has sometimes problems with masses, it give carbon mass of calcium
                # this code makes sure that the masses are correct:
                def strip_numbers_from_atom_name(atom_name):
                    return re.match("([a-zA-Z]+)", atom_name).group(0)

                for atom in range(len(syst_fingerprint_heavy_atoms.atoms)):
                    if np.abs(syst_fingerprint_heavy_atoms.atoms[atom].mass - reordered_system.atoms[atom].mass) > 1:
                        if (strip_numbers_from_atom_name(
                                syst_fingerprint_heavy_atoms.atoms[atom].name) == strip_numbers_from_atom_name(
                            reordered_system.atoms[atom].name)):
                            syst_fingerprint_heavy_atoms.atoms[atom].mass = reordered_system.atoms[atom].mass

                if np.linalg.norm(syst_fingerprint_heavy_atoms.positions - reordered_system.positions)>0.01:
                    rmsd_fp = rmsd(syst_fingerprint_heavy_atoms.positions, reordered_system.positions, superposition = True)
                else:
                    rmsd_fp = 0.0

                if rmsd_fp < best_rmsd:
                    logger.info(f"Found better permutation: {rmsd_fp:} < {best_rmsd:}")
                    best_rmsd = rmsd_fp
                    best_mapping = reversed_concatenated_mapping

                if best_rmsd < 0.01:  # if rmsd_fp is small there is no point of searching further
                    break
 
    upper_names = np.array([name.upper() for name in syst_fingerprint_heavy_atoms.atoms.names])
    temp = np.where([upper_names == metal_name.upper()])[0]

    # Reconstruct full mapping, including hydrogen bonds

    # we create full maps of all atoms:
    # G_fingerprint = nx.Graph(
    #    MDAnalysis.topology.guessers.guess_bonds(syst_fingerprint.atoms[1:], syst_fingerprint.atoms[1:].positions,
    #                                             vdwradii={metal_name: 1.0}))  # TODO assumption that first atom is the metal

    G_fingerprint = syst_to_graph(syst_fingerprint_pbc.atoms[1:], vdwradii={metal_name: 1.0})

    # no_metal = connected_cut_system_pbc.select_atoms(f"not index {metal_index:d}", sorted=False)
    no_metal = connected_cut_system_pbc[1:]  # .select_atoms(f"not index {metal_index:d}", sorted=False)
    G_site = syst_to_graph(no_metal.atoms, vdwradii={metal_name: 1.0})
    # G_site = nx.Graph(
    #    MDAnalysis.topology.guessers.guess_bonds(no_metal.atoms, no_metal.atoms.positions))

    # we copy indexes of heavy atoms
    heavy = {node: node for node in G_fingerprint.nodes if node in G_fingerprint_heavy_atoms.nodes}
    # to all atoms not present (i.e., hydrogens) we add -1:
    hydrogens = {node: -1 for node in G_fingerprint.nodes if node not in G_fingerprint_heavy_atoms.nodes}

    # the atoms which are not part of mapping are hydrogens:
    inverted_best_mapping = dict((b, a) for a, b in best_mapping.items())
    not_mapped = {atom.index: -1 for atom in no_metal.atoms if atom.index not in inverted_best_mapping}

    # assigning the indexes to the graphs:
    nx.set_node_attributes(G_fingerprint, {**heavy, **hydrogens}, "idx")
    nx.set_node_attributes(G_site, {**inverted_best_mapping, **not_mapped}, "idx")

    iso = isomorphism.GraphMatcher(G_fingerprint, G_site, node_match=lambda n1, n2: n1['idx'] == n2['idx'])

    if iso.is_isomorphic():
        # we copy one of the mappings, now we do not have to worry about hydrogen atoms:
        whole_best_mapping = iso.mapping
    else:
        print("Error, not found the mapping")
        raise

    if len(temp) == 1:
        whole_best_mapping[temp[0]] = metal_index
    else:
        logger.info("ERROR, more than one metal in the site")
        raise

    # Remove the end atoms from the mapping:
    logger.info(f"Best mapping:")
    # logger.info(f"    Mapping: {best_mapping}")
    logger.info(f"    Mapping: {whole_best_mapping}")

    mapping_end_atoms = []
    for index1 in whole_best_mapping:
        # logger.info(
        #    f"      Map: {index1} -> {best_mapping[index1]} ({syst_fingerprint.atoms[index1]} -> {cage.atoms[best_mapping[index1]]})")
        if whole_best_mapping[index1] in end_atoms:
            mapping_end_atoms.append(index1)
            logger.info("          End atom! Will be removed")

    logger.info(f"    RMSD: {best_rmsd:}")
    logger.info(f"    Removing end atoms: {end_atoms:} in fingerprint nonation: {mapping_end_atoms:}")
    for end_atom in mapping_end_atoms:
        whole_best_mapping.pop(end_atom)



    return whole_best_mapping, best_rmsd

    '''
    mapping_end_atoms = []
    for index1 in best_mapping:
        #logger.info(
        #    f"      Map: {index1} -> {best_mapping[index1]} ({syst_fingerprint.atoms[index1]} -> {cage.atoms[best_mapping[index1]]})")
        if best_mapping[index1] in end_atoms:
            mapping_end_atoms.append(index1)
            logger.info("          End atom! Will be removed")

    logger.info(f"    RMSD: {best_rmsd:}")
    logger.info(f"    Removing end atoms: {end_atoms:} in fingerprint nonation: {mapping_end_atoms:}")
    for end_atom in mapping_end_atoms:
        best_mapping.pop(end_atom)

    return best_mapping, best_rmsd
    '''
