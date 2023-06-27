import argparse
import MDAnalysis
import networkx as nx
import re
import numpy as np

from MDAnalysis.lib.distances import distance_array
from copy import deepcopy

import rdkit
from rdkit import Chem
from MDAnalysis.topology.guessers import guess_types

from networkx import isomorphism
import parmed as pmd


try:
    from cgbind2pmd.log import logger
    from cgbind2pmd.mapping import map_two_structures, unwrap
    from cgbind2pmd.utils import new_directory, mdanalysis_to_rdkit,strip_numbers_from_atom_names
except:
    from log import logger
    from mapping import map_two_structures, unwrap
    from utils import new_directory, mdanalysis_to_rdkit, strip_numbers_from_atom_names


def find_metal_indices(cage, metal_name):
    if not hasattr(cage.atoms[0], 'element'):
        guessed_elements = MDAnalysis.topology.guessers.guess_types(strip_numbers_from_atom_names(cage.atoms.names))
        cage.universe.add_TopologyAttr('elements', guessed_elements)
    metal_indices = [idx for idx, element in enumerate(cage.atoms.elements) if
                     element.title() == metal_name.title()]
    n_metals = len(metal_indices)
    
    if n_metals == 0: # this happens when guesser fails, for example Pd becomes P. We need to try using just names 
        metal_indices = [idx for idx, name in enumerate(cage.atoms.names) if
                         name[:len(metal_name)].title() == metal_name.title()]
        n_metals = len(metal_indices)    
    
    return metal_indices, n_metals

# TODO find_bound_whole_ligands and find_bound_ligands_nx are the same, aren't they?
# TODO THey are the same, and remove this after sevela weeks, if nothing breaks (19/05/2023)

'''
def find_bound_whole_ligands(metal, cage, G_sub_ligands, cutoff_covalent=3):    
    Find ligands which are bound to the metal

    metal: mdanalysis
    

    bound_ligands = []
    closest_atoms = []

    for ligand in G_sub_ligands:
        closest_atoms_ligand = []

        metal_ligand_distances = distance_array(cage[list(ligand.nodes)].positions, metal.position)

        # arg_min = np.argmin(metal_ligand_distances)

        for arg_min in range(len(metal_ligand_distances)):

            if metal_ligand_distances[arg_min] < cutoff_covalent:
                closest_atoms_ligand.append(cage[list(ligand.nodes)][arg_min].index)

                if ligand not in bound_ligands:
                    bound_ligands.append(ligand)
        if len(closest_atoms_ligand) > 0:
            closest_atoms.append(closest_atoms_ligand)

    n_ligands = len(bound_ligands)
    return bound_ligands, closest_atoms

'''

def find_bound_ligands_nx(cage, metal_index, cutoff=7, cutoff_covalent=3.0, closest_neighbhors = 3, neighbhor_cutoff=None):
    '''
    Finds bound ligands to metal, assumes that atoms within cutoff_covalent (default 3) are bound to metal.
    Returns list of sub-graphs of bound ligands; it does not cut ligands (so they can be uneven)

    :param metal_index:
    :param cutoff:
    :param cutoff_covalent:
    :return:
    '''

    if cutoff is not None:
        cut_sphere = cage.select_atoms(f'around {cutoff:f} index {metal_index:d}')
    else:
        cut_sphere = cage.select_atoms(f"not index {metal_index:d}")

    additional_atom_types = {}
    for atom_type in list(set(cut_sphere.types)):
        if atom_type not in MDAnalysis.topology.tables.vdwradii:
            additional_atom_types[atom_type]=0.1

    metal = cage.atoms[metal_index]

    G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(cut_sphere.atoms, cut_sphere.atoms.positions, box=cage.dimensions, vdwradii=additional_atom_types))
    nx.set_node_attributes(G_cage, {atom.index: atom.name[0] for atom in cut_sphere.atoms}, "name")
    G_sub_cages = [G_cage.subgraph(a) for a in nx.connected_components(G_cage)]

    closest_atoms_string = f"Closest atoms around {metal_index:d}:"
    closest_atoms_ligands = []

    G_sub_cages_bound = []
    for G_sub_cage in G_sub_cages:
        closest_atoms = []  # closest atoms, we assume that they are donors of electrons
        clusters_of_atoms = cage.atoms[list(G_sub_cage)]


        all_metal_cluster_distances = distance_array(clusters_of_atoms.positions, metal.position, box=cage.dimensions)
        if np.min(all_metal_cluster_distances) < cutoff_covalent:
            G_sub_cages_bound.append(G_sub_cage)
            G_indices = list(G_sub_cage.nodes)


            # we now search closest atoms, if atoms are futher then closest_neighbhors then we also add it to clostest atoms
            # this is needed for bidente lignads
            ordered = np.argsort(all_metal_cluster_distances.T[0])
            close_atoms = []
            close_atoms.append(ordered[0])

            for order in ordered[1:]:
                if all_metal_cluster_distances.T[0][order] < cutoff_covalent: #order has to be change to real number!

                    append = True
                    for temp_atom in close_atoms:
                        if nx.shortest_path_length(G_sub_cage, G_indices[temp_atom], G_indices[order]) < closest_neighbhors:
                            append = False
                    if append == True:
                        close_atoms.append(order)

            #temp_atom = clusters_of_atoms.atoms[np.argmin(all_metal_cluster_distances)]
            for atom in close_atoms:
                temp_atom = clusters_of_atoms.atoms[atom]
                closest_atoms_string+=f" {temp_atom.name:s} {temp_atom.index:d}: {np.min(all_metal_cluster_distances):f} A"
                closest_atoms.append(temp_atom.index)
                
        
            closest_atoms_ligands.append(closest_atoms)

    logger.info(closest_atoms_string)

    G_sub_cages_bound_cutoff = []
    closest_atoms_ligands_cutoff = []
    if neighbhor_cutoff is not None:
        for idx, closest_atoms in enumerate(closest_atoms_ligands):
            nodes = []
            for closest_atom in closest_atoms:
                # we take radius from the closest atom (so not the metal!) that is why it is -1
                new_nodes = nx.generators.ego_graph(G_sub_cages_bound[idx], closest_atom, radius=neighbhor_cutoff - 1).nodes
                nodes += list(new_nodes)
            nodes = list(set(nodes))

            for G_sub_cage in nx.connected_components(nx.subgraph(G_sub_cages_bound[idx], nodes)):
                closest_atoms_ligand = []
                
                for close_atom in closest_atoms:
                    if close_atom in list(G_sub_cage):
                        closest_atoms_ligand.append(close_atom)

                G_sub_cages_bound_cutoff.append(G_sub_cages_bound[idx].subgraph(G_sub_cage))
                closest_atoms_ligands_cutoff.append(closest_atoms_ligand)

        G_sub_cages_bound = G_sub_cages_bound_cutoff
        closest_atoms_ligands = closest_atoms_ligands_cutoff


    #ourput is list with n_ligands graphs representing diffrent ligands, and closest atoms (donors of electrons)
    return G_sub_cages_bound, closest_atoms_ligands

def find_closest_and_add_rings(metal_index, cage, bound_ligands, selected_closest_atoms, aromaticity):

    #  # , cutoff=7, cutoff_covalent=3.0, closest_neighbhors=3)
    # TODO THIS NEEDs to be checked for the knots/cages

    site_link_atoms_single = []
    site_single = []
    site_ligands_single =[]

    for idx_bound_ligand in range(len(bound_ligands)):
        close_atoms_and_in_rings = []
        selected_bound_ligand = nx.Graph(bound_ligands[idx_bound_ligand])
        ##selected_closest_atom = selected_closest_atoms[a]

        # the old way, to be remove0d # TODO (19/05/2023)
        # cut_sphere = cage[[metal_index] + list(selected_bound_ligand.nodes)].select_atoms(
        #    f'around {cutoff:f} index {metal_index:d}')

        # we select neighbours of 3, which are always added
        neighbour_indices = []
        # we also want to add rings directly connected to donor atoms (so not only rings which consist of a donor atom)
        direct_neighbour_indices = []
        for closest_atom in selected_closest_atoms[idx_bound_ligand]:
            if closest_atom in selected_bound_ligand.nodes:
                neighbour_indices += list(nx.generators.ego_graph(selected_bound_ligand, closest_atom, radius=2).nodes)
                direct_neighbour_indices += list(
                    nx.generators.ego_graph(selected_bound_ligand, closest_atom, radius=1).nodes)

        cut_sphere = cage[[metal_index] + neighbour_indices]
        cut_sphere_direct = cage[[metal_index] + direct_neighbour_indices]

        atoms_within_cutoff = list(cut_sphere.indices)
        atoms_close_to_donor = list(cut_sphere_direct.indices)

        # Find rings with the colvalently bound atom:
        rings_with_closest_atom = []

        rings = nx.cycle_basis(selected_bound_ligand)

        aromatic_rings = []

        def recursive_rings(ring_of_intrest, rings, aromaticity):
            # Search recursively for the largest connected aromatic system (we extend one by one the system)
            for ring in rings:
                # Check if rings have intersection and assert that ring is not only a subset of the bigger ring
                if len(np.intersect1d(ring_of_intrest, ring)) >= 2 and len(set(ring_of_intrest + ring)) > len(
                        ring_of_intrest) and len(set(ring_of_intrest + ring)) > len(ring):
                    new_ring = list(set(ring + ring_of_intrest))
                    if all(aromaticity[new_ring]):
                        return recursive_rings(new_ring, rings, aromaticity)
            return ring_of_intrest

        for idx1, ring in enumerate(rings):
            # we extend ring to the largest connected aromatic system
            new_ring = recursive_rings(ring, rings, aromaticity)
            # since small rings can be part of the same larger aromatic system we check if not already added:
            if new_ring not in aromatic_rings:
                aromatic_rings.append(new_ring)

        for ring in aromatic_rings:
            # for selected_closest_atom in selected_closest_atoms[idx_bound_ligand]:
            for selected_closest_atom in atoms_close_to_donor:
                if selected_closest_atom in ring:
                    rings_with_closest_atom += ring
                else:
                    # We remove other rings, which would become fragmets otherwise
                    for atom in ring:
                        if atom in atoms_within_cutoff:
                            atoms_within_cutoff.remove(atom)

        # We merge two criteria
        close_atoms_and_in_rings = list(set(atoms_within_cutoff + rings_with_closest_atom))

        # we remove atoms which are alone, not connected to anything
        # this might be effect that ring was removed, but its hydrogen left behind, as a result, the fragment of the ring is reconstructured, which we don't want
        Gsub = selected_bound_ligand.subgraph(close_atoms_and_in_rings)
        close_atoms_and_in_rings = list(np.concatenate(
            [list(Gsub_connect) for Gsub_connect in nx.connected_components(Gsub) if len(Gsub_connect) > 1]))
        close_atoms_plus_adj = close_atoms_and_in_rings[:]

        adj_atoms = []

        # we add adjecent atoms
        for edge in selected_bound_ligand.edges:
            if edge[0] in close_atoms_and_in_rings and edge[1] not in close_atoms_and_in_rings:
                adj_atoms.append(edge[1])
            elif edge[1] in close_atoms_and_in_rings and edge[0] not in close_atoms_and_in_rings:
                adj_atoms.append(edge[0])

        close_atoms_plus_adj += adj_atoms

        # This does not work when for example things are charged, like [O-] at the end, remove in fuuture (26/05/2023) TODO
        # Link atoms are all atoms which are not hydrogens and which are added:
        # link_atoms = [atom.index for atom in cage[adj_atoms] if atom.type != 'H']
        # site_link_atoms += link_atoms

        # This is how it should be done We search for link atoms, link atoms will have different node degree in
        # subgraph of select site then whole ligand
        Gsub_extended = selected_bound_ligand.subgraph(close_atoms_plus_adj)
        degrees_subset = nx.degree(Gsub_extended)
        degrees_full = nx.degree(selected_bound_ligand)

        site_link_atoms_single += [node for node in dict(degrees_subset) if degrees_subset[node] != degrees_full[node]]

        site_single += close_atoms_plus_adj

        site_ligands_single.append(close_atoms_plus_adj)

        '''
        site_link_atoms += [node for node in dict(degrees_subset) if degrees_subset[node] != degrees_full[node]]
        site += close_atoms_plus_adj
        site_ligands.append(close_atoms_plus_adj)
        '''

    return site_single, site_ligands_single, site_link_atoms_single


def find_closest_and_add_rings_iterate(metal_indices, cage, all_metal_indecies):
    '''
    Select closest atoms connected to metal.


    :param metal_indices:       metal_indecies
    :param cage:                MDAnalysis input
    :param G_sub_ligands:
    :return:
    '''
    # This is the main function, we cut the thing, and add the rings
    cutoff = 4
    unsaturated_atoms = []
    # site = [metal_index]
    # site = [0]
    # site_ligands = []

    binding_sites_graphs = []
    binding_sites_site_link_atoms= []
    

    # Use MDAnalysis (interface to RDKIT) to get aromaticity
    for metal_index in all_metal_indecies:
        cage.atoms[metal_index].type = 'H' # we change type of metal for hydrogen, becasue vdwradii is required for aromaticity


    if not hasattr(cage.atoms[0], 'element'):
        guessed_elements = MDAnalysis.topology.guessers.guess_types(strip_numbers_from_atom_names(cage.names))
        for metal_index in metal_indices:
            guessed_elements[metal_index] = 'H'
        cage.universe.add_TopologyAttr('elements', guessed_elements)

    aromaticity = MDAnalysis.topology.guessers.guess_aromaticities(cage)


    for metal_index in metal_indices:

        site = []
        #site_ligands = []
        site_link_atoms = []

        #metal = cage.atoms[metal_index]
        #bound_ligands2, selected_closest_atoms2 = find_bound_whole_ligands(metal, cage, G_sub_ligands)

        bound_ligands, selected_closest_atoms = find_bound_ligands_nx(cage, metal_index,cutoff=None)
        site_single, _ , site_link_atoms_single = find_closest_and_add_rings(metal_index, cage,bound_ligands, selected_closest_atoms, aromaticity)
        site+= site_single
        #site_ligands += site_ligands_single
        site_link_atoms += site_link_atoms_single

        # Add metal to the selection, metal goes first
        binding_sites_graphs.append([metal_index] + list(sorted(set(site))))  # The sorting is required for the later checking of uniqueness (making ordering of atoms in figerprint hassle-free)
        binding_sites_site_link_atoms.append(site_link_atoms)

        # find unsaturated atomsmap_two_structures
        # we check if the degree of the node is the same for cutted ligand and the whole ligand
        # if diffrent, than this atoms is unsaturated
        '''
        subgraph = selected_bound_ligand.subgraph(close_atoms_plus_adj)

        for degree in subgraph.degree():
            if degree[1] != selected_bound_ligand.degree()[degree[0]]:
                unsaturated_atoms.append(degree[0])
        '''


    return binding_sites_graphs, binding_sites_site_link_atoms


def check_uniqueness(binding_sites_graphs, site_link_atoms, structure, metal_name, rmsd_cutoff=2):
    # we check if the sites exist in the file

    unique_sites = []
    unique_site_link_atoms = []

    logger.info(f"Uniqueness check: Number of binding sites to check: {len(binding_sites_graphs):}")
    print(f"Uniqueness check: Number of binding sites to check: {len(binding_sites_graphs):}")
    

    for site, site_link_atom in zip(binding_sites_graphs, site_link_atoms):
        exist = False
        print("next site, already unique sites:", len(unique_sites))

        for site_fingerprint in unique_sites:
            logger.info(f"\tChecking uniqueness of site and fingerprint with atoms:{len(site):}, {len(site_fingerprint):}")
            print(f"\tChecking uniqueness of site:{len(site):}, {len(site_fingerprint):}")

            # the structures are sorted because MDAnalysis selection atoms sorts them
            #_, rmsd = map_two_structures(site[0], structure[np.sort(site)], structure[np.sort(site_fingerprint)], metal_name)

            _, rmsd = map_two_structures(0, structure[site], structure[site_fingerprint],
                                         metal_name)




            logger.info(f"RMSD between current structure and already check structures: {rmsd:}")
            print(f"RMSD between current structure and already check structures: {rmsd:}")
            if rmsd < rmsd_cutoff:
                print("yep!")
                exist = True
                break

        if exist == False:
            print("Adding unique site")
            unique_sites.append(site)
            unique_site_link_atoms.append(site_link_atom)

    logger.info(f"Number of unique sites: {len(unique_sites):}")
    return unique_sites, unique_site_link_atoms


from tempfile import mkdtemp
import openbabel

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def add_hydrogens(selection_site, metal_name, add_atoms_to_this_atom):
    '''
    Adds hydrogen to unsaturated carbon

    :param selection_site:
    :param metal_name:
    :param add_atoms_to_this_atom:
    :return:
    '''
    logger.info("Adding hydrogens")
    # we just align the methane and add the 3 hydrogens from it, optimises with PBE0/def2-SVP
    methane = np.array([[-0.00000935761277, -0.00000001492783, 0.00002482854051],
                        [-0.00004478038360, 0.00000000997987, 1.09805688107415],
                        [1.03524344997147, 0.00000000862120, -0.36607166084943],
                        [-0.51759466202540, -0.89654704617818, -0.36600502948382],
                        [-0.51759464994970, 0.89654704250494, -0.36600501928140]])

    amonia = np.array([[ 0.        ,  0.        ,  0.        ],
                       [ 1.01734118,  0.06528425,  0.02511132],
                       [-0.30613666,  0.79279584, -0.56356069],
                       [-0.30612386,  0.21132302,  0.94946643]])

    water = np.array([[ 0.        ,  0.        ,  0.        ],
                      [ 0.9611851 , -0.02275931, -0.04795393],
                      [-0.27035532, -0.39612253, -0.83468673]])

    sulfide = np.array([[ 0.        ,  0.        ,  0.        ],
                        [ 1.33344014, -0.14807463, -0.13973411],
                        [-0.25250322, -0.96369492, -0.90939079]])

    positions = []
    for link_atom in add_atoms_to_this_atom:

        if selection_site.universe.dimensions is not None:
            G_selection_site = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(selection_site.atoms,
                                                                                 selection_site.atoms.positions,
                                                                                 vdwradii={metal_name: 3},
                                                                                 box = selection_site.universe.dimensions))
        else:
            G_selection_site = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(selection_site.atoms,
                                                                                 selection_site.atoms.positions,
                                                                                 vdwradii={metal_name: 3}))

        vector = list(G_selection_site.edges(link_atom))

        if len(vector) > 1:
            print("Error, more than one atom connected")
            raise

        atom_connected = vector[0][0]
        if atom_connected == link_atom:
            atom_connected = vector[0][1]

        atom1_pos = selection_site.select_atoms(f'index {link_atom:d}')[0].position
        atom2_pos = selection_site.select_atoms(f'index {atom_connected:d}')[0].position

        bond_vector = atom2_pos - atom1_pos
        bond_vector /= np.linalg.norm(vector)

        molecule = methane # we guess it is carbon, unless it is something from N, O, S, other elements not included

        link_atom_name = selection_site.select_atoms(f'index {link_atom:d}')[0].name[0]
        if link_atom_name == 'N':
            molecule = amonia
        elif link_atom_name == 'O':
            molecule = water
        elif link_atom_name == 'S':
            molecule = sulfide

        rot = rotation_matrix_from_vectors(molecule[1] - molecule[0], bond_vector)
        newpos = rot.dot(molecule.T).T[2:]

        positions.extend(np.array(newpos) + atom1_pos)

    # create the Universe
    n_atoms = len(positions)
    logger.info(f"Number of atoms after adding the hydrogens {n_atoms:}")

    if n_atoms > 0:
        hydrogens = MDAnalysis.Universe.empty(n_atoms, trajectory=True)
        hydrogens.add_TopologyAttr('name', ['H'] * n_atoms)
        hydrogens.add_TopologyAttr('type', ['H'] * n_atoms)
        hydrogens.atoms.positions = positions
        new_syst = MDAnalysis.Merge(selection_site.atoms, hydrogens.atoms)
    else:
        new_syst = MDAnalysis.Merge(selection_site.atoms) #this makes sub-universe

    renumered_add_atoms_to_this_atom = [idx for idx, name in enumerate(selection_site.atoms.indices) if name in add_atoms_to_this_atom]

    return new_syst, renumered_add_atoms_to_this_atom

def find_atom_to_ligand_membership(new_syst, metal_name):
    new_site_no_metal = new_syst.select_atoms(f"not name {metal_name.title():s} {metal_name.upper():s} {metal_name.lower():s}")
    
    if new_site_no_metal.universe.dimensions is not None:
        G_all_short_ligands = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(new_site_no_metal.atoms, new_site_no_metal.atoms.positions, box=new_site_no_metal.universe.dimensions))
    else:
        G_all_short_ligands = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(new_site_no_metal.atoms, new_site_no_metal.atoms.positions))

    nx.set_node_attributes(G_all_short_ligands, {atom.index: atom.name[0] for atom in new_site_no_metal.atoms},
                           "name")
    G_sub_short_ligands = [G_all_short_ligands.subgraph(a) for a in nx.connected_components(G_all_short_ligands)]

    # Here we just take the subgraphs, so selected atoms, but nx makes the order random, this is not a problem general
    # unless you want to take the extra hydrogens out
    ligands_atoms_membership = [list(G_sub_short_ligands[a].nodes) for a in range(len(G_sub_short_ligands))]

    return ligands_atoms_membership

def find_ligand_pattern(new_syst, ligands_nodes):
    '''
    Check if the site is composed of the same ligand fragments (we want to symmetrize paramters for the same sites)

    :param new_syst:
    :param ligands_nodes:
    :return:
    '''

    #new_syst2 = new_syst.copy()


    unique_ligands = []
    unique_ligands_pattern = []

    guessed_elements = guess_types(strip_numbers_from_atom_names(new_syst.atoms.names))
    new_syst.add_TopologyAttr('elements', guessed_elements)

    for ligands_node in ligands_nodes:
        exist = False

        selected_ligand_1 = MDAnalysis.Merge(new_syst.atoms[ligands_node])
        #selected_ligand_1.atoms.write("temp1.pdb")
        #mol1 = Chem.MolFromPDBFile("temp1.pdb", removeHs=False) # remove (?) 13/06 TODO

        mol1 = selected_ligand_1.atoms.convert_to("RDKIT") # TODO remove
        #mol1 = mdanalysis_to_rdkit(selected_ligand_1)

        #if mol1 != mol1b: # TODO remove
        #    raise


        for idx, unique_lingad in enumerate(unique_ligands):

            selected_ligand_2 = MDAnalysis.Merge(new_syst.atoms[unique_lingad])
            #selected_ligand_2.atoms.write("temp2.pdb")# remove (?) 13/06 TODO
            #mol2 = Chem.MolFromPDBFile("temp2.pdb", removeHs=False)

            mol2 = selected_ligand_2.atoms.convert_to("RDKIT") # TODO remove
            #mol2 = mdanalysis_to_rdkit(selected_ligand_2)

            #if mol2 !=mol2b:
            #    raise




            if (mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1)):
                exist = True
                unique_ligands_pattern.append(idx)

                break

        if exist == False:
            unique_ligands_pattern.append(len(unique_ligands))
            unique_ligands.append(ligands_node)

    return unique_ligands_pattern, unique_ligands

def renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern, extra_atoms, link_atoms):
    '''
    We renumber the atoms that there are togheter. Also, identical ligands also have to have the same numbering

    :param new_syst:
    :param metal_name:
    :return:
    '''

    new_link_atoms = []
    new_extra_atoms = []

    metals = new_syst.select_atoms(f"name {metal_name.title():s} {metal_name.upper():s} {metal_name.lower():s}")
    for metal in metals:
        metal.element = metal_name.title()


    new_ligands = [metals]
    
    n_total_atoms = len(metals)

    for idx, ligands_node in enumerate(ligands_atoms_membership):
        print(idx)
        sorted_extra_atoms = [idx for idx, a in enumerate(ligands_node) if a in extra_atoms]
        sorted_link_atoms = [idx for idx, a in enumerate(ligands_node) if a in link_atoms]

        #new_syst.atoms[ligands_node].write("temp1.pdb") # TODO check why do we need this, remove (?) 13/06 TODO
        #new_syst.atoms[ligands_node].write("temp1.xyz") # xyz has better precision then pdb, important for the strain calculations. We might even try to omit that TODO
        #selected_ligand_1 = MDAnalysis.Universe("temp1.xyz") # atom indices are fixed and cannot be change, this is a simple way of renumbering themvremove (?) 13/06 TODO
        #selected_ligand_1b = new_syst.atoms[ligands_node] #remove (?) 13/06 TODOremove (?) 13/06 TODO
        selected_ligand_1 = MDAnalysis.Merge(new_syst.atoms[ligands_node]) # Merge is needed to renumber things
        #selected_ligand_1 = new_syst.atoms[ligands_node]
        #mol1 = mdanalysis_to_rdkit(selected_ligand_1)
        #print([atom.GetChiralTag() for atom in mol1.GetAtoms()])

        #mol1 = Chem.MolFromPDBFile("temp1.pdb", removeHs=False)
        mol1 = selected_ligand_1.atoms.convert_to("RDKIT") # TODO I probably should use openbabel to guess eleemtns

        G_ligand_1 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(selected_ligand_1.atoms, selected_ligand_1.atoms.positions))
        nx.set_node_attributes(G_ligand_1, {atom.index: atom.name[0] for atom in selected_ligand_1.atoms}, "name")
        nx.set_node_attributes(G_ligand_1, {atom.GetChiralTag() for atom in mol1.GetAtoms()}, "chirality")

        #mol1b = mdanalysis_to_rdkit(selected_ligand_1) # TODO remove
        #if mol1 != mol1b:
        #    raise

        #new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]].write("temp2.pdb") remove (?) 13/06 TODO
        #new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]].write("temp2.xyz")remove (?) 13/06 TODO
        #selected_ligand_2 = MDAnalysis.Universe("temp2.pdb")remove (?) 13/06 TODO
        #selected_ligand_2 = MDAnalysis.Universe("temp2.xyz")remove (?) 13/06 TODO

        #selected_ligand_2 = new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]]
        selected_ligand_2 = MDAnalysis.Merge(new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]])

        #mol1 = mdanalysis_to_rdkit(selected_ligand_2)

        #mol1 = Chem.MolFromPDBFile("temp2.pdb", removeHs=False)
        mol1 = selected_ligand_2.atoms.convert_to("RDKIT")
        G_ligand_2 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(selected_ligand_2.atoms, selected_ligand_2.atoms.positions))
        nx.set_node_attributes(G_ligand_2, {atom.index: atom.name[0] for atom in selected_ligand_2.atoms}, "name")
        nx.set_node_attributes(G_ligand_2, {atom.GetChiralTag() for atom in mol1.GetAtoms()}, "chirality")

        print([atom.GetChiralTag() for atom in mol1.GetAtoms()])

        new_ligand = selected_ligand_2.copy()


        #mol1b = mdanalysis_to_rdkit(selected_ligand_2)
        #if mol1 != mol1b:
        #    raise


        iso = isomorphism.GraphMatcher(G_ligand_1, G_ligand_2,
                                       node_match=lambda n1, n2: n1['name'] == n2['name'] and n1['chirality'] == n2[
                                           'chirality'])
        if iso.is_isomorphic():
            for node in iso.mapping:
                new_ligand.atoms[iso.mapping[node]].position = selected_ligand_1.atoms[node].position

            link_atoms_ligand = []
            for link_atom in sorted_link_atoms:
                #link_atoms.append(n_total_atoms + iso.mapping[link_atom])
                link_atoms_ligand.append(iso.mapping[link_atom])
            new_link_atoms.append(link_atoms_ligand)

            extra_atoms_ligand = []
            for extra_atom in sorted_extra_atoms:
                #new_extra_atoms.append(n_total_atoms + iso.mapping[extra_atom])
                extra_atoms_ligand.append(iso.mapping[extra_atom])
            new_extra_atoms.append(extra_atoms_ligand)

            print(new_ligands)
            new_ligands.append(new_ligand.atoms)

            n_total_atoms += len(new_ligand.atoms)
            print(new_ligands, n_total_atoms)
        else:
            logger.info(f"The graphs are not isomorphic, index: {idx:}")

    return new_ligands, new_extra_atoms, new_link_atoms


def extract_metal_structure(filename, topol_filename, metal_name, output=None, check_uniquness=True, all_metal_names=None):
    '''
    It takes a structure and tries to find structures around metals.

    :param filename:
    :param metal_name:
    :param output:
    :param check_uniquness: it is used for the strain calculations
    additional_metals -> metals which are present in topology but removed for parametrization of small molecules
    :return:
    '''
    syst = MDAnalysis.Universe(filename)
    topol = pmd.load_file(topol_filename)

    cage = syst.atoms

    metal_indices, n_metals = find_metal_indices(syst, metal_name)

    all_metals_indices = []
    if all_metal_names is not None:
        for additional_metal in all_metal_names:
            all_metals_indices += find_metal_indices(syst, additional_metal)[0]
    else:
        all_metals_indices = metal_indices

    if n_metals == 0:
        logger.info(f"Metal {metal_name:} not found")
        raise

    logger.info(f"Metal indices to check {metal_indices:} metal name:{metal_name:}")
    all_ligands_atoms = cage.select_atoms(f"not index {' '.join(map(str, all_metals_indices)):}")

    G_all_ligands = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(all_ligands_atoms.atoms, all_ligands_atoms.atoms.positions, box=all_ligands_atoms.dimensions))
    nx.set_node_attributes(G_all_ligands, {atom.index: atom.name[0] for atom in all_ligands_atoms.atoms}, "name")

    G_sub_ligands = [G_all_ligands.subgraph(a) for a in nx.connected_components(G_all_ligands)]

    binding_sites_graphs, site_link_atoms = find_closest_and_add_rings_iterate(metal_indices, cage, all_metals_indices)

    if check_uniquness:
        unique_sites, unique_site_link_atoms = check_uniqueness(binding_sites_graphs, site_link_atoms,  cage, metal_name)
    else:
        unique_sites = binding_sites_graphs
        unique_site_link_atoms = site_link_atoms

    metal_sites = []



    for n_site, unique_site in enumerate(unique_sites):
        new_directory(f"{output:s}{n_site:d}")

        site_topol = deepcopy(topol)
        strip_atoms = [idx for idx in list(cage.indices) if idx not in unique_site]
        site_topol.strip(f"@{','.join(list(map(str, np.array(strip_atoms) + 1))):s}")



        selection_site = cage[unique_site]
        #cage.atoms.write("cage.pdb") # TODO not sure why this is here
        #selection_site.atoms.write(f"temp_temp{n_site:}.pdb") #TODO remove (?)

        if selection_site.universe.dimensions is not None:
            selection_site.atoms[0].type = metal_name # previously, we changed it to H (aromaticity) we change it back
            selection_site = unwrap(selection_site, metal_name)

        new_syst, new_link_atoms = add_hydrogens(selection_site, metal_name, unique_site_link_atoms[n_site])

        #Hydrogens are added at the end of the file
        extra_atoms = new_syst.atoms[len(selection_site):].atoms.indices
        #new_syst.atoms.write("temp.gro")

        # we split site, and we group atoms for diffrent ligands
        ligands_atoms_membership = find_atom_to_ligand_membership(new_syst, metal_name)

        # we check if some ligands are the same
        unique_ligands_pattern, unique_ligands = find_ligand_pattern(new_syst, ligands_atoms_membership)

        # save the ligands, guess it's charges using RDKIT
        unique_ligands_charges = []
        unique_ligand_filenames = []
        unique_ligands_smiles = []





        unique_ligand_topols = []

        for idx, unique_ligand in enumerate(unique_ligands):
            ligand_topol = deepcopy(site_topol)

            #strip_atoms = [idx for idx in list(new_syst.atoms[:len(selection_site)].indices) if idx not in unique_ligand]

            selected_atoms_indices = [idx for idx in new_syst.atoms[:len(selection_site)].indices if idx in unique_ligand]
            ligand_topol.strip(f"!@{','.join(list(map(str, np.array(selected_atoms_indices) + 1))):s}")
            unique_ligand_topols.append(ligand_topol)

            ligand_coord = MDAnalysis.Merge(new_syst.atoms[unique_ligand])
            ligand_coord.atoms.write(f"{output:s}{n_site:d}/ligand_{idx:}.pdb") # TODO do we need both (?)
            ligand_coord.atoms.write(f"{output:s}{n_site:d}/ligand_{idx:}.xyz")
            
            unique_ligand_filenames.append(f"{output:s}{n_site:d}/ligand_{idx:}.xyz")

            # Guess charges of the ligand fragment, whatch out, this usually works, but sometimes gives wierd result!
            mol = ligand_coord.atoms.convert_to("RDKIT")

            #mol = MDAnalysis.Merge(new_syst.atoms[ligands_node])

            #mol1b = mdanalysis_to_rdkit(ligand_coord) # TODO remove (?)
            #if mol != mol1b:
            #    raise

            #mol = mdanalysis_to_rdkit(ligand_coord)

            charge = rdkit.Chem.GetFormalCharge(mol)

            unique_ligands_charges.append(charge)
            unique_ligands_smiles.append(rdkit.Chem.MolToSmiles(mol))

        # TODO rename? sorted here means categorised to the specific ligand
        renumbered_ligands, sorted_extra_atoms, sorted_link_atoms = renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern, extra_atoms, new_link_atoms)
        
        # Parmed needs that the same residues are grouped togheter:
        sorted_renumbered_ligands = [renumbered_ligands[a] for a in np.append(0, np.argsort(unique_ligands_pattern) + 1)]
        #sorted_renumbered_extra_atoms = [sorted_extra_atoms[a] for a in np.argsort(unique_ligands_pattern))]
        #sorted_renumbered_extra_link_atoms = [sorted_link_atoms[a] for a in np.append(0, np.argsort(unique_ligands_pattern))]


        #sorted_extra_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if a in extra_atoms]
        #sorted_link_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if a in new_link_atoms]

        #sorted_link_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if
        #                      a in unique_site_link_atoms]

        # Saving the files ------------------------
        new_cage = MDAnalysis.Merge(*sorted_renumbered_ligands)  # , dimensions=crystal.dimensions)
        new_cage.atoms.write(f"{output:s}{n_site:d}/site.pdb")
        new_cage.atoms.write(f"{output:s}{n_site:d}/site.xyz")

        metal_topol = deepcopy(site_topol)
        # metal is the first (we renumbered them in selection site), in amber counting starts from 1
        metal_topol.strip("!@1")
        new_site_topol = metal_topol
        for unique in unique_ligands_pattern:
            # renumbering of the atoms is done according to the first ligand so it should be exactly as for the coordination
            new_site_topol+=unique_ligand_topols[unique]


        # MDAnalysis writer captials (i.g., FE), orca needs title type (i.g., Fe)
        File = open(f"{output:s}{n_site:d}/site.xyz")
        text = File.read()
        File.close()
        
        with open(f"{output:s}{n_site:d}/site.xyz", "w") as File:
            File.write(text.title())

        starting_index = [0]
        for struct in sorted_renumbered_ligands:
            starting_index.append(starting_index[-1] + len(struct))

        indecies = [list(range(starting_index[a - 1], starting_index[a])) for a in range(1, len(starting_index))]

        charge_pattern = [unique_ligands_charges[unique_ligands_pattern[a]] for a in np.argsort(unique_ligands_pattern)]
        smiles_pattern = [unique_ligands_smiles[unique_ligands_pattern[a]] for a in np.argsort(unique_ligands_pattern)]
        sorted_renumbered_extra_atoms = [np.array(sorted_extra_atoms[a]) for a in np.argsort(unique_ligands_pattern)]
        sorted_renumbered_extra_atoms = list(np.concatenate([starting_index[idx + 1] + extra_atoms for idx, extra_atoms in enumerate(sorted_renumbered_extra_atoms)]))
        
        sorted_renumbered_extra_link_atoms = [np.array(sorted_link_atoms[a]) for a in np.argsort(unique_ligands_pattern)]
        sorted_renumbered_extra_link_atoms = list(np.concatenate([starting_index[idx + 1] + link_atoms for idx, link_atoms in
                                         enumerate(sorted_renumbered_extra_link_atoms)]))

        with open("INFO.dat", "a") as File:
            File.write(f"ligand_pattern:{','.join(list(map(str,np.sort(unique_ligands_pattern)))):}\n")
            File.write(f"link_atoms:{','.join(list(map(str, sorted_renumbered_extra_link_atoms))):}\n")
            File.write(f"extra_atoms:{','.join(list(map(str, sorted_renumbered_extra_atoms))):}\n")
            File.write(f"starting_index:{','.join(list(map(str, starting_index))):}\n")
            File.write(f"charge_pattern:{','.join(list(map(str, charge_pattern))):}\n")
            File.write(f"smiles_pattern:{','.join(list(map(str, smiles_pattern))):}\n")

        selected_metal_site_list = [metal_name, 0, 1, f"{output:s}{n_site:d}", unique_ligand_filenames,
                                    np.sort(unique_ligands_pattern), sorted_renumbered_extra_link_atoms,
                                    sorted_renumbered_extra_atoms, starting_index, indecies, charge_pattern,
                                    smiles_pattern, new_site_topol]
        # charge is set up to redicouls value, to make sure that I remember to change it, #TODO chage it to 0

        metal_sites.append(selected_metal_site_list)
    
    return metal_sites
    
    


def read_info_file(self, filename='INFO.dat'):  # TODO this should be everywhere, od subrutines should be removed
    File = open(filename)
    text = File.read()
    File.close()

    n_sites = text.count("ligand_pattern")
    print("Number of sites", n_sites)

    unique_ligands_pattern = []
    link_atoms = []
    additional_atoms = []
    starting_indecies = []
    charge_pattern = []
    indecies = []

    for line in text.splitlines():
        if "ligand_pattern:" in line:
            unique_ligands_pattern.append(list(map(int, line[15:].split(','))))

        elif "link_atoms:" in line:
            if ',' in line:
                link_atoms.append(list(map(int, line[11:].split(','))))
            else:
                link_atoms.append([])

        elif "extra_atoms:" in line:
            if ',' in line:
                additional_atoms.append(list(map(int, line[12:].split(','))))
            else:
                additional_atoms.append([])

        elif "starting_index:" in line:
            starting_index = list(map(int, line[15:].split(',')))
            starting_indecies.append(starting_index)
            indecies.append(
                [list(range(starting_index[a - 1], starting_index[a])) for a in range(1, len(starting_index))])

        elif "charge_pattern:" in line:
            charge_pattern.append(list(map(int, line[15:].split(','))))

    sites = []
    for n_site in range(n_sites):
        filename = f"site{n_site:d}.xyz"

        sites.append(metal_site(metal_name, metal_charge, filename, unique_ligands_pattern[n_site], link_atoms[n_site],
                                additional_atoms[n_site], starting_indecies[n_site], indecies[n_site],
                                charge_pattern[n_site]))
    self.unique_sites = sites


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-o", help="Clean structure")
    parser.add_argument("-metal_name", help="Name of the metal")
    return parser.parse_args()



if __name__ == '__main__':
    args = get_args()
    if args.o is not None:
        extract_metal_structure(args.f, args.metal_name, args.o)


