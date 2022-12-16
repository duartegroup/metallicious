import os

import networkx as nx
import MDAnalysis
from itertools import permutations
from MDAnalysis.analysis import align
from networkx.algorithms import isomorphism
import numpy as np
import re

import MDAnalysis.transformations

try:
    from cgbind2pmd.log import logger
except:
    from log import logger



def map_two_structures(metal_index, connected_cut_system, syst_fingerprint, metal_name, end_atoms=[]):
    print(os.getcwd())
    G_fingerprint = nx.Graph(
        MDAnalysis.topology.guessers.guess_bonds(syst_fingerprint.atoms, syst_fingerprint.atoms.positions, vdwradii={metal_name:1.0}))
    nx.set_node_attributes(G_fingerprint, {atom.index: atom.name[0] for atom in syst_fingerprint.atoms}, "name")
    G_fingerprint_subs = [G_fingerprint.subgraph(a) for a in nx.connected_components(G_fingerprint)]

    # Can we cut here?
    # do we need this (?)
    # cut_sphere = self.cage.select_atoms(f'index {metal_index:d} or around {cutoff:f} index {metal_index:d}')
    # connected_cut_system = cut_sphere.select_atoms("index " + " ".join(map(str, selected_atoms)))
    # I think this should work just fine, no? TODO:


    #first solve pbd problem
    new_universe = connected_cut_system.universe.copy()
    G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(connected_cut_system.atoms,
                                                               connected_cut_system.atoms.positions,
                                                               vdwradii={connected_cut_system[0].type: 3},
                                                               box=new_universe.dimensions))
    bonds = list(G_cage.edges)
    new_universe.add_TopologyAttr('bonds', bonds)
    connected_cut_system_pbc = new_universe.atoms[connected_cut_system.atoms.indices]
    transform = MDAnalysis.transformations.unwrap(connected_cut_system_pbc)
    new_universe.trajectory.add_transformations(transform)
    # solved


    no_metal = connected_cut_system_pbc.select_atoms(f"not index {metal_index:d}")
    G_site = nx.Graph(
        MDAnalysis.topology.guessers.guess_bonds(no_metal.atoms, no_metal.atoms.positions))
    nx.set_node_attributes(G_site, {atom.index: atom.name[0] for atom in no_metal.atoms}, "name")
    G_site_subs = [G_site.subgraph(a) for a in nx.connected_components(G_site)]

    number_ligands_bound = len(G_site_subs)

    best_rmsd = 1e10
    best_mapping = None

    indecies2number = {}
    for a, b in enumerate(connected_cut_system_pbc.indices):
        indecies2number[b] = a

    permutations_list = permutations(range(number_ligands_bound))

    # We check all permutation of residues
    for perm in list(permutations_list):
        mappings = []
        # some metal sites are not symmetric, if they are not we don't check them
        length_match = True

        # we check all permutation of renumering of atoms in resudes
        for a, b in enumerate(perm):
            # if lenght of ligands is diffrent, than this is deffinitly wrong permutation
            if len(G_fingerprint_subs[a]) != len(G_site_subs[b]):
                mappings = []
                length_match = False
                break

            iso = isomorphism.GraphMatcher(G_fingerprint_subs[a], G_site_subs[b],
                                           node_match=lambda n1, n2: n1['name'] == n2['name'])
            mappings.append([subgraph_mapping for subgraph_mapping in iso.subgraph_isomorphisms_iter()])

        if length_match:
            # we tabule all the prmutations of permutations (renumering atoms)
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

            for c in range(len(all_possible_mappings)):
                new_dict = {}
                for ligand in all_possible_mappings[c]:
                    new_dict.update(ligand)

                # Reorder, left hand side
                new_new_dict = {}
                for d in sorted(new_dict.keys()):
                    new_new_dict[d] = new_dict[d]

                # TODO We assume that metal is first

                reordered_system = connected_cut_system_pbc.atoms[
                    [0] + [indecies2number[value] for value in new_new_dict.values()]]  # Find where metal is
                #make sure that the type is preseved (esptially challening with names starting with C and N)
                reordered_system.atoms[0].type = metal_name
                reordered_system.atoms[0].mass = 0.0

                '''
                #Solvig PBC problem:
                new_universe = reordered_system.universe.copy()
                G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(reordered_system.atoms, reordered_system.atoms.positions,
                                                                           vdwradii={reordered_system[0].type: 3},
                                                                           box=new_universe.dimensions))
                bonds = list(G_cage.edges)
                new_universe.add_TopologyAttr('bonds', bonds)
                reordered_system_pbc = new_universe.atoms[reordered_system.atoms.indices]
                transform = MDAnalysis.transformations.unwrap(reordered_system_pbc)
                new_universe.trajectory.add_transformations(transform)
                # solved
                '''

                # MDAnalysis has sometimes problems with masses, it give carbon mass of calcium
                # this code makes sure that the masses are correct:
                def strip_numbers_from_atom_name(atom_name):
                    return re.match("([a-zA-Z]+)", atom_name).group(0)

                for atom in range(len(syst_fingerprint.atoms)):
                    if np.abs(syst_fingerprint.atoms[atom].mass - reordered_system.atoms[atom].mass) > 1:
                        if (strip_numbers_from_atom_name(
                                syst_fingerprint.atoms[atom].name) == strip_numbers_from_atom_name(
                                reordered_system.atoms[atom].name)):
                            syst_fingerprint.atoms[atom].mass = reordered_system.atoms[atom].mass



                _, rmsd = align.alignto(syst_fingerprint, reordered_system)

                if rmsd > 10: #if RMSD is large probably is due to the PBC problem... this takes much longer then
                    new_universe = reordered_system.universe.copy()
                    G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(reordered_system.atoms,
                                                                               reordered_system.atoms.positions,
                                                                               vdwradii={reordered_system[0].type: 3},
                                                                               box=new_universe.dimensions))
                    bonds = list(G_cage.edges)
                    new_universe.add_TopologyAttr('bonds', bonds)
                    reordered_system_pbc = new_universe.atoms[reordered_system.atoms.indices]
                    transform = MDAnalysis.transformations.unwrap(reordered_system_pbc)
                    new_universe.trajectory.add_transformations(transform)
                    _, rmsd2 = align.alignto(syst_fingerprint, reordered_system)
                    # solved


                if rmsd < best_rmsd:
                    print(rmsd, "<", best_rmsd)
                    best_rmsd = rmsd
                    best_mapping = new_new_dict

                if best_rmsd < 0.01: #if rmsd is really small there is no point of searching further
                    break

    upper_names = np.array([name.upper() for name in syst_fingerprint.atoms.names])
    temp = np.where([upper_names == metal_name.upper()])[0]
    if len(temp) == 1:
        best_mapping[temp[0]] = metal_index
    else:
        logger.info("ERROR, more than one metal in the site")
        raise

    # Remove the end atoms from the mapping:

    logger.info(f"Best mapping:")
    logger.info(f"    Mapping: {best_mapping}")

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
