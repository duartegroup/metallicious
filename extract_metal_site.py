import argparse
import MDAnalysis
import networkx as nx
import re
import numpy as np

from MDAnalysis.lib.distances import distance_array

from main import map_two_structures

import rdkit
from rdkit import Chem
from MDAnalysis.topology.guessers import guess_types

from networkx import isomorphism


def find_metal_indices(cage, metal_name):
    metal_indices = [a for a, name in enumerate(cage.atoms.names) if
                     name[:len(metal_name)].title() == metal_name.title()]
    n_metals = len(metal_indices)
    return metal_indices, n_metals


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


def find_closest_and_add_rings(metal_indices, cage, G_sub_ligands ):
    # This is the main function, we cut the thing, and add the rings
    cutoff = 4
    unsaturated_atoms = []
    # site = [metal_index]
    # site = [0]
    # site_ligands = []

    binding_sites_graphs = []

    for metal_index in metal_indices:

        site = [metal_index]
        site_ligands = []
        site_link_atoms = []

        metal = cage.atoms[metal_index]
        bound_ligands, selected_closest_atoms = find_bound_whole_ligands(metal,cage, G_sub_ligands)

        for a in range(len(bound_ligands)):
            close_atoms_and_in_rings = []
            selected_bound_ligand = nx.Graph(bound_ligands[a])
            ##selected_closest_atom = selected_closest_atoms[a]

            cut_sphere = cage[[metal_index] + list(selected_bound_ligand.nodes)].select_atoms(
                f'around {cutoff:f} index {metal_index:d}')
            atoms_within_cutoff = list(cut_sphere.indices)

            # Find rings with the colvalently bound atom:
            rings_with_closest_atom = []

            for ring in nx.cycle_basis(selected_bound_ligand):
                for selected_closest_atom in selected_closest_atoms[a]:
                    if selected_closest_atom in ring:
                        rings_with_closest_atom += ring
                    else:
                        # We remove other rings, which would become fragmets otherwise
                        for atom in ring:
                            if atom in atoms_within_cutoff:
                                atoms_within_cutoff.remove(atom)

            # We merge two crateria
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

        binding_sites_graphs.append(list(set(site)))

        # find unsaturated atoms
        # we check if the degree of the node is the same for cutted ligand and the whole ligand
        # if diffrent, than this atoms is unsaturated
        '''
        subgraph = selected_bound_ligand.subgraph(close_atoms_plus_adj)

        for degree in subgraph.degree():
            if degree[1] != selected_bound_ligand.degree()[degree[0]]:
                unsaturated_atoms.append(degree[0])
        '''
        return binding_sites_graphs, site_link_atoms


def check_uniqueness(binding_sites_graphs, structure, metal_name):
    # we check if the sites exist in the file

    unique_sites = []


    for site in binding_sites_graphs:
        exist = False

        for site_fingerprint in unique_sites:
            print("Checking", len(site), len(site_fingerprint))

            _, rmsd = map_two_structures(0, structure[np.sort(site)], structure[np.sort(site_fingerprint)], metal_name)
            print("RMSD", rmsd)
            if rmsd < 1.0:
                exist = True
                break

        if exist == False:
            unique_sites.append(site)

    print("we have this amoutnt of unique sites:", len(unique_sites))
    return unique_sites


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
    print(n_atoms)
    if n_atoms > 0:
        hydrogens = MDAnalysis.Universe.empty(n_atoms, trajectory=True)
        hydrogens.add_TopologyAttr('name', ['H'] * n_atoms)
        hydrogens.add_TopologyAttr('type', ['H'] * n_atoms)
        hydrogens.atoms.positions = positions
        new_syst = MDAnalysis.Merge(selection_site.atoms, hydrogens.atoms)
    else:
        new_syst = MDAnalysis.Merge(selection_site.atoms) #this makes sub-universe

    return new_syst

def find_atom_to_ligand_membership(new_syst, metal_name):
    new_site_no_metal = new_syst.select_atoms(f"not name {metal_name.title():s} {metal_name.upper():s} {metal_name.lower():s}")
    if new_site_no_metal.universe.dimensions is not None:
        print("A")
        G_all_short_ligands = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(new_site_no_metal.atoms, new_site_no_metal.atoms.positions, box=new_site_no_metal.universe.dimensions))
    else:
        print("B")
        G_all_short_ligands = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(new_site_no_metal.atoms, new_site_no_metal.atoms.positions))

    nx.set_node_attributes(G_all_short_ligands, {atom.index: atom.name[0] for atom in new_site_no_metal.atoms},
                           "name")
    G_sub_short_ligands = [G_all_short_ligands.subgraph(a) for a in nx.connected_components(G_all_short_ligands)]
    ligands_atoms_membership = [list(G_sub_short_ligands[a].nodes) for a in range(len(G_sub_short_ligands))]
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

        for a, unique_lingad in enumerate(unique_ligands):

            selected_ligand_2 = new_syst.atoms[unique_lingad]
            selected_ligand_2.atoms.write("temp.pdb")
            mol2 = Chem.MolFromPDBFile("temp.pdb")

            if (mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1)):
                exist = True
                unique_ligands_pattern.append(a)

                break

        if exist == False:
            unique_ligands_pattern.append(len(unique_ligands))
            unique_ligands.append(ligands_node)

    return unique_ligands_pattern, unique_ligands

def renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern):
    '''
    We renumber the atoms that there are togheter. Also, identical ligands also have to have the same numbering

    :param new_syst:
    :param metal_name:
    :return:
    '''

    metals = new_syst.select_atoms(f"name {metal_name.title():s} {metal_name.upper():s} {metal_name.lower():s}")
    for metal in metals:
        metal.element = metal_name.title()
    new_ligands = [metals]

    for idx, ligands_node in enumerate(ligands_atoms_membership):
        new_syst.atoms[ligands_node].write("temp1.pdb")
        selected_ligand_1 = MDAnalysis.Universe("temp1.pdb")
        mol1 = Chem.MolFromPDBFile("temp1.pdb")
        G_ligand_1 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(selected_ligand_1.atoms, selected_ligand_1.atoms.positions))
        nx.set_node_attributes(G_ligand_1, {atom.index: atom.name[0] for atom in selected_ligand_1.atoms}, "name")
        nx.set_node_attributes(G_ligand_1, {atom.GetChiralTag() for atom in mol1.GetAtoms()}, "chirality")

        new_syst.atoms[unique_ligands[unique_ligands_pattern[idx]]].write("temp2.pdb")
        selected_ligand_2 = MDAnalysis.Universe("temp2.pdb")
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
            new_ligands.append(new_ligand.atoms)
        else:
            print(idx, "The graphs are not isomorphic")

    return new_ligands


def extract_metal_structure(filename,metal_name, output=None):
    '''
    The main loop, it takes a structure and tries to find structures around metals.

    :param filename:
    :param metal_name:
    :param output:
    :return:
    '''
    syst = MDAnalysis.Universe(filename)


    cage = syst.atoms

    metal_indices, n_metals = find_metal_indices(syst, metal_name)

    print(metal_indices, metal_name)

    all_ligands_atoms = cage.select_atoms("not index " + " ".join(map(str,metal_indices)))

    G_all_ligands = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(all_ligands_atoms.atoms, all_ligands_atoms.atoms.positions, box=all_ligands_atoms.dimensions))
    nx.set_node_attributes(G_all_ligands, {atom.index: atom.name[0] for atom in all_ligands_atoms.atoms}, "name")

    G_sub_ligands = [G_all_ligands.subgraph(a) for a in nx.connected_components(G_all_ligands)]

    binding_sites_graphs, site_link_atoms = find_closest_and_add_rings(metal_indices, cage, G_sub_ligands)

    unique_sites = check_uniqueness(binding_sites_graphs, cage, metal_name)
    print("Site link", site_link_atoms)


    for idx_site, unique_site in enumerate(unique_sites):

        selection_site = cage[unique_site]
        print("unique_site", unique_site)
        cage.atoms.write("cage.pdb")
        selection_site.atoms.write(f"temp_temp{idx_site:}.pdb")

        # wrapping around pbc #TODO separate function (?)
        if selection_site.universe.dimensions is not None:
            print("dim=",selection_site.universe.dimensions)

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


        new_syst = add_hydrogens(selection_site, metal_name, site_link_atoms)


        #Hydrogens are added at the end of the file
        extra_atoms = new_syst.atoms[len(selection_site):].atoms.indices
        #new_syst.atoms.write("temp.gro")

        # we split site, and we group atoms for diffrent ligands
        ligands_atoms_membership = find_atom_to_ligand_membership(new_syst, metal_name)

        # we check if some of the ligands are the same
        unique_ligands_pattern, unique_ligands = find_ligand_pattern(new_syst, ligands_atoms_membership)

        # save the ligands
        for idx, unique_ligand in enumerate(unique_ligands):
            new_syst.atoms[unique_ligand].write(f"ligand{idx_site:d}_{idx:}.pdb")
            new_syst.atoms[unique_ligand].write(f"ligand{idx_site:d}_{idx:}.xyz")

        renumbered_ligands = renumer_ligands(new_syst, metal_name, ligands_atoms_membership, unique_ligands, unique_ligands_pattern)

        # Parmed needs that the same residues are grouped togheter:
        sorted_renumbered_ligands = [renumbered_ligands[a] for a in np.append(0, np.argsort(unique_ligands_pattern) + 1)]
        sorted_extra_atoms = [idx for idx, a in enumerate(
            np.concatenate([ligands_atoms_membership[a] for a in np.argsort(unique_ligands_pattern)])) if a in extra_atoms]


        # Saving the files ------------------------
        new_cage = MDAnalysis.Merge(*sorted_renumbered_ligands)  # , dimensions=crystal.dimensions)
        new_cage.atoms.write(f"{output:s}{idx_site:d}.pdb")
        new_cage.atoms.write(f"{output:s}{idx_site:d}.xyz")

        # MDAnalysis writer captials (i.g., FE), orca needs title type (i.g., Fe)
        File = open(f"{output:s}{idx_site:d}.xyz")
        text = File.read()
        File.close()
        with open(f"{output:s}{idx_site:d}.xyz", "w") as File:
            File.write(text.title())

        starting_index = [0]
        for struct in sorted_renumbered_ligands:
            starting_index.append(starting_index[-1] + len(struct))

        with open("INFO.dat", "a") as File:
            File.write(f"ligand_pattern:{','.join(list(map(str,np.sort(unique_ligands_pattern)))):}\n")
            File.write(f"link_atoms:{','.join(list(map(str, sorted_extra_atoms))):}\n")
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


