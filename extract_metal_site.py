import argparse
import MDAnalysis
import networkx as nx
import re
import numpy as np

from MDAnalysis.lib.distances import distance_array


import rdkit
from rdkit import Chem
from MDAnalysis.topology.guessers import guess_types

from networkx import isomorphism


try:
    from cgbind2pmd.log import logger
    from cgbind2pmd.mapping import map_two_structures

except:
    from log import logger
    from mapping import map_two_structures


def find_metal_indices(cage, metal_name):
    metal_indices = [a for a, name in enumerate(cage.atoms.names) if
                     name[:len(metal_name)].title() == metal_name.title()]
    n_metals = len(metal_indices)
    return metal_indices, n_metals

# TODO find_bound_whole_ligands and find_bound_ligands_nx are the same, aren't they?
# TODO THey are the same, and remove this after sevela weeks, if nothing breaks (19/05/2023)
'''
def find_bound_whole_ligands(metal, cage, G_sub_ligands, cutoff_covalent=3):
    '''
    Find ligands which are bound to the metal

    metal: mdanalysis
    '''

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

def find_bound_ligands_nx(cage, metal_index, cutoff=7, cutoff_covalent=3.0, closest_neighbhors = 3):
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
        cut_sphere = cage
    #metal_type = cut_sphere.select_atoms(f'index {metal_index:}').atoms[0].type


    # TODO this does not work for the connected ligands (??)
    G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(cut_sphere.atoms, cut_sphere.atoms.positions, box=cage.dimensions))#, vdwradii={metal_type:cutoff_covalent}))
    nx.set_node_attributes(G_cage, {atom.index: atom.name[0] for atom in cut_sphere.atoms}, "name")
    G_sub_cages = [G_cage.subgraph(a) for a in nx.connected_components(G_cage)]

    closest_atoms_string = f"Closest atoms around {metal_index:d}:"
    closest_atoms_ligands = []

    G_sub_cages_bound = []
    for G_sub_cage in G_sub_cages:
        closest_atoms = []  # closest atoms, we assume that they are donors of electrons
        clusters_of_atoms = cage.atoms[list(G_sub_cage)]
        metal = cage.atoms[metal_index]

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

    #ourput is list with n_ligands graphs representing diffrent ligands, and closest atoms (donors of electrons)
    return G_sub_cages_bound, closest_atoms_ligands




def find_closest_and_add_rings(metal_indices, cage, G_sub_ligands):
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
    for metal_index in metal_indices:
        cage.atoms[metal_index].type = 'H' # we change type of metal for hydrogen, becasue vdwradii is required for aromaticity

    if not hasattr(cage.atoms[0], 'element'):
        guessed_elements = MDAnalysis.topology.guessers.guess_types(cage.names)
        for metal_index in metal_indices:
            guessed_elements[metal_index] = 'H'
        cage.universe.add_TopologyAttr('elements', guessed_elements)

    aromaticity = MDAnalysis.topology.guessers.guess_aromaticities(cage)


    for metal_index in metal_indices:

        site = []
        site_ligands = []
        site_link_atoms = []

        metal = cage.atoms[metal_index]
        bound_ligands2, selected_closest_atoms2 = find_bound_whole_ligands(metal, cage, G_sub_ligands)
        bound_ligands, selected_closest_atoms = find_bound_ligands_nx(cage, metal_index, cutoff=None)#, cutoff=7, cutoff_covalent=3.0, closest_neighbhors=3)
        # TODO THIS NEEDs to be checked for the knots/cages

        for idx_bound_ligand in range(len(bound_ligands)):
            close_atoms_and_in_rings = []
            selected_bound_ligand = nx.Graph(bound_ligands[idx_bound_ligand])
            ##selected_closest_atom = selected_closest_atoms[a]


            # the old way, to be removed # TODO
            #cut_sphere = cage[[metal_index] + list(selected_bound_ligand.nodes)].select_atoms(
            #    f'around {cutoff:f} index {metal_index:d}')

            # we select neighbours of 3
            neighbour_indices = []
            for closest_atom in selected_closest_atoms[idx_bound_ligand]:
                if closest_atom in selected_bound_ligand.nodes:
                    neighbour_indices += list(nx.generators.ego_graph(selected_bound_ligand, closest_atom, radius=2).nodes)

            cut_sphere = cage[[metal_index] + neighbour_indices]


            atoms_within_cutoff = list(cut_sphere.indices)

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
                for selected_closest_atom in selected_closest_atoms[idx_bound_ligand]:
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
            # (this might be effect that ring was removed, but its hydrogen left behind, as a result, the fragment of the ring is reconstructured, which we don't want)
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

            close_atoms_plus_adj+=adj_atoms

            #Link atoms are all atoms which are not hydrogens and which are added:
            link_atoms = [atom.index for atom in cage[adj_atoms] if atom.type != 'H']
            site_link_atoms += link_atoms


            site += close_atoms_plus_adj
            site_ligands.append(close_atoms_plus_adj)

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


def check_uniqueness(binding_sites_graphs, site_link_atoms, structure, metal_name):
    # we check if the sites exist in the file

    unique_sites = []
    unique_site_link_atoms = []

    logger.info(f"Uniqueness check: Number of binding sites to check: {len(binding_sites_graphs):}")
    for site, site_link_atom in zip(binding_sites_graphs, site_link_atoms):
        exist = False

        for site_fingerprint in unique_sites:
            logger.info("\tChecking uniqueness of site:", len(site), len(site_fingerprint))

            # the structures are sorted because MDAnalysis selection atoms sorts them
            #_, rmsd = map_two_structures(site[0], structure[np.sort(site)], structure[np.sort(site_fingerprint)], metal_name)

            _, rmsd = map_two_structures(0, structure[site], structure[site_fingerprint],
                                         metal_name)

            logger.info("RMSD between current structure and already check structures:", rmsd)
            if rmsd < 1.0:
                exist = True
                break

        if exist == False:
            unique_sites.append(site)
            unique_site_link_atoms.append(site_link_atom)

    logger.info("Number of unique sites:", len(unique_sites))
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

        rot = rotation_matrix_from_vectors(methane[1] - methane[0], bond_vector)
        newpos = rot.dot(methane.T).T[2:]

        positions.extend(np.array(newpos) + atom1_pos)

    # create the Universe
    n_atoms = len(positions)
    logger.info("Number of atoms after adding the hydrogens", n_atoms)

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
    ligands_atoms_membership = [list(G_sub_short_ligands[a].nodes) for a in range(len(G_sub_short_ligands))] # TODO I made it sorted, check if this helped! well. it helps in the ordering... it breakes others
    #ligands_atoms_membership = [sorted(list(G_sub_short_ligands[a].nodes)) for a in range(len(G_sub_short_ligands))] # TODO I made it sorted, check if this helped! well. it helps in the ordering... it breakes others
    return ligands_atoms_membership

def find_ligand_pattern(new_syst, ligands_nodes):
    '''
    Check if the site is composed of the same ligand fragments (we want to symmetrize paramters for the same sites)

    :param new_syst:
    :param ligands_nodes:
    :return:
    '''


    unique_ligands = []
    unique_ligands_pattern = []

    guessed_elements = guess_types(new_syst.atoms.names)
    new_syst.add_TopologyAttr('elements', guessed_elements)

    for ligands_node in ligands_nodes:
        exist = False

        selected_ligand_1 = new_syst.atoms[ligands_node]
        selected_ligand_1.atoms.write("temp.pdb")
        mol1 = Chem.MolFromPDBFile("temp.pdb")

        for idx, unique_lingad in enumerate(unique_ligands):

            selected_ligand_2 = new_syst.atoms[unique_lingad]
            selected_ligand_2.atoms.write("temp.pdb")
            mol2 = Chem.MolFromPDBFile("temp.pdb")

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

        sorted_extra_atoms = [idx for idx, a in enumerate(ligands_node) if a in extra_atoms]
        sorted_link_atoms = [idx for idx, a in enumerate(ligands_node) if a in link_atoms]

        new_syst.atoms[ligands_node].write("temp1.pdb")
        new_syst.atoms[ligands_node].write("temp1.xyz") # xyz has better precision then pdb, important for the strain calculations. We might even try to omit that TODO
        selected_ligand_1 = MDAnalysis.Universe("temp1.xyz")

        mol1 = Chem.MolFromPDBFile("temp1.pdb")
        G_ligand_1 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(selected_ligand_1.atoms, selected_ligand_1.atoms.positions))
        nx.set_node_attributes(G_ligand_1, {atom.index: atom.name[0] for atom in selected_ligand_1.atoms}, "name")
        nx.set_node_attributes(G_ligand_1, {atom.GetChiralTag() for atom in mol1.GetAtoms()}, "chirality")

        new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]].write("temp2.pdb")
        new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]].write("temp2.xyz")
        #selected_ligand_2 = MDAnalysis.Universe("temp2.pdb")
        selected_ligand_2 = MDAnalysis.Universe("temp2.xyz")
        mol1 = Chem.MolFromPDBFile("temp2.pdb")
        G_ligand_2 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(selected_ligand_2.atoms, selected_ligand_2.atoms.positions))
        nx.set_node_attributes(G_ligand_2, {atom.index: atom.name[0] for atom in selected_ligand_2.atoms}, "name")
        nx.set_node_attributes(G_ligand_2, {atom.GetChiralTag() for atom in mol1.GetAtoms()}, "chirality")

        new_ligand = selected_ligand_2.copy()

        iso = isomorphism.GraphMatcher(G_ligand_1, G_ligand_2,
                                       node_match=lambda n1, n2: n1['name'] == n2['name'] and n1['chirality'] == n2[
                                           'chirality'])
        if iso.is_isomorphic():
            for node in iso.mapping:
                new_ligand.atoms[iso.mapping[node]].position = selected_ligand_1.atoms[node].position
            for link_atom in sorted_link_atoms:
                new_link_atoms.append(n_total_atoms + iso.mapping[link_atom])
            for extra_atom in sorted_extra_atoms:
                new_extra_atoms.append(n_total_atoms + iso.mapping[extra_atom])

            new_ligands.append(new_ligand.atoms)
            n_total_atoms += len(new_ligand.atoms)
        else:
            logger.info(f"The graphs are not isomorphic, index: {idx:}")

    return new_ligands, new_extra_atoms, new_link_atoms


def extract_metal_structure(filename, metal_name, output=None, check_uniquness=True):
    '''
    It takes a structure and tries to find structures around metals.

    :param filename:
    :param metal_name:
    :param output:
    :param check_uniquness: it is used for the strain calculations
    :return:
    '''
    syst = MDAnalysis.Universe(filename)

    cage = syst.atoms

    metal_indices, n_metals = find_metal_indices(syst, metal_name)

    logger.info("Metal indices to check {metal_indices:} metal name:{metal_name:}")

    all_ligands_atoms = cage.select_atoms("not index " + " ".join(map(str,metal_indices)))

    G_all_ligands = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(all_ligands_atoms.atoms, all_ligands_atoms.atoms.positions, box=all_ligands_atoms.dimensions))
    nx.set_node_attributes(G_all_ligands, {atom.index: atom.name[0] for atom in all_ligands_atoms.atoms}, "name")

    G_sub_ligands = [G_all_ligands.subgraph(a) for a in nx.connected_components(G_all_ligands)]

    binding_sites_graphs, site_link_atoms = find_closest_and_add_rings(metal_indices, cage, G_sub_ligands)

    if check_uniquness:
        unique_sites, unique_site_link_atoms = check_uniqueness(binding_sites_graphs, site_link_atoms,  cage, metal_name)
    else:
        unique_sites = binding_sites_graphs
        unique_site_link_atoms = site_link_atoms


    for n_site, unique_site in enumerate(unique_sites):

        selection_site = cage[unique_site]
        cage.atoms.write("cage.pdb")
        selection_site.atoms.write(f"temp_temp{n_site:}.pdb")

        # wrapping around pbc #TODO separate function (?)
        if selection_site.universe.dimensions is not None:
            logger.info(f"Dimensionality {selection_site.universe.dimensions:}")

            G_selection_site = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(selection_site.atoms,
                                                                                 selection_site.atoms.positions,
                                                                                 vdwradii={metal_name: 3},
                                                                                 box=selection_site.universe.dimensions))
            bonds = list(G_selection_site.edges)
            selection_site.universe.add_TopologyAttr('bonds', bonds)
            #connected_cut_system_pbc = selection_site.universe.atoms[connected_cut_system.atoms.indices]
            transform = MDAnalysis.transformations.unwrap(selection_site)
            selection_site.universe.trajectory.add_transformations(transform)
            # solved


        new_syst, new_link_atoms = add_hydrogens(selection_site, metal_name, unique_site_link_atoms[n_site])


        #Hydrogens are added at the end of the file
        extra_atoms = new_syst.atoms[len(selection_site):].atoms.indices
        #new_syst.atoms.write("temp.gro")

        # we split site, and we group atoms for diffrent ligands
        ligands_atoms_membership = find_atom_to_ligand_membership(new_syst, metal_name)

        # we check if some of the ligands are the same
        unique_ligands_pattern, unique_ligands = find_ligand_pattern(new_syst, ligands_atoms_membership)

        # save the ligands
        for idx, unique_ligand in enumerate(unique_ligands):
            new_syst.atoms[unique_ligand].write(f"ligand{n_site:d}_{idx:}.pdb")
            new_syst.atoms[unique_ligand].write(f"ligand{n_site:d}_{idx:}.xyz")

        renumbered_ligands, sorted_extra_atoms, sorted_link_atoms  = renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern, extra_atoms, new_link_atoms)

        # Parmed needs that the same residues are grouped togheter:
        sorted_renumbered_ligands = [renumbered_ligands[a] for a in np.append(0, np.argsort(unique_ligands_pattern) + 1)]
        #sorted_extra_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if a in extra_atoms]
        #sorted_link_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if a in new_link_atoms]

        #sorted_link_atoms = [idx for idx, a in enumerate(
        #    np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if
        #                      a in unique_site_link_atoms]

        # Saving the files ------------------------
        new_cage = MDAnalysis.Merge(*sorted_renumbered_ligands)  # , dimensions=crystal.dimensions)
        new_cage.atoms.write(f"{output:s}{n_site:d}.pdb")
        new_cage.atoms.write(f"{output:s}{n_site:d}.xyz")

        # MDAnalysis writer captials (i.g., FE), orca needs title type (i.g., Fe)
        File = open(f"{output:s}{n_site:d}.xyz")
        text = File.read()
        File.close()
        with open(f"{output:s}{n_site:d}.xyz", "w") as File:
            File.write(text.title())

        starting_index = [0]
        for struct in sorted_renumbered_ligands:
            starting_index.append(starting_index[-1] + len(struct))

        with open("INFO.dat", "a") as File:
            File.write(f"ligand_pattern:{','.join(list(map(str,np.sort(unique_ligands_pattern)))):}\n")
            File.write(f"link_atoms:{','.join(list(map(str, sorted_link_atoms))):}\n")
            File.write(f"extra_atoms:{','.join(list(map(str, sorted_extra_atoms))):}\n")
            File.write(f"starting_index:{','.join(list(map(str, starting_index))):}\n")



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


