import logging

import MDAnalysis
import networkx as nx
import numpy as np

from networkx.algorithms import isomorphism
import parmed as pmd
from copy import deepcopy

from metallicious.antechamber_interface import antechamber, check_antechamber_if_available
from metallicious.extract_metal_site import find_metal_indices
from metallicious.initial_site import create_metal_topol
from metallicious.mapping import unwrap
from metallicious.log import logger


def mapping_itp_coords(syst, ligand_file):
    if ligand_file is None:
        return False, False, False, False

    G_lig = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(syst.atoms, syst.atoms.positions))
    nx.set_node_attributes(G_lig, {atom.index: atom.name[0] for atom in syst.atoms}, "name")
    topol = pmd.load_file(ligand_file)
    bonds = [(bond.atom1.idx, bond.atom2.idx) for bond in topol.bonds]
    G_top = nx.Graph(bonds)
    nx.set_node_attributes(G_top, {atom.idx: atom.name[0] for atom in topol.atoms}, "name")
    iso = nx.isomorphism.GraphMatcher(G_lig, G_top, node_match=lambda n1, n2: n1['name'] == n2['name'])
    if not iso.is_isomorphic():
        return False, False, False, False
    else:
        ordered = MDAnalysis.Universe.empty(len(topol.atoms), trajectory=True)
        ordered.add_TopologyAttr('name')
        for a, _ in enumerate(ordered.atoms):
            ordered.atoms[a].name = topol.atoms[a].name

        return True, topol, ordered, G_top

def prepare_initial_topology(filename, metal_names, metal_charge, output_coord, output_top, metal_vdw, ligand_topol=None):
    if check_antechamber_if_available()== False:
        return ImportError("Not antechamber detected (preparation of initial topology is done by antechamber)")



    crystal = MDAnalysis.Universe(filename)
    crystal.residues.resids = 0 # Antechamber requires one residue

    table = MDAnalysis.topology.tables.vdwradii
    table['Cl'] = table['CL']
    metal_ions = None

    if len(metal_names) > 0:
        if crystal.dimensions is not None:
            for metal_name in metal_names:
                crystal = unwrap(crystal, metal_name, metal_cutoff=0.1)

        #metal_ions = crystal.select_atoms(f"name {metal_name.title():s}* {metal_name.upper():s}* {metal_name.lower():s}* ")
        metal_indecies = []
        for metal_name in metal_names:
            indecies, _ = find_metal_indices(crystal, metal_name)
            metal_indecies.append(indecies)

        metal_ions = crystal.atoms[np.concatenate(metal_indecies)]
        ligands = crystal.atoms - metal_ions

        #ligands = crystal.select_atoms(f"not name {metal_name.title():s}* {metal_name.upper():s}* {metal_name.lower():s}* ")
        G1 = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table,
                                                     box=crystal.dimensions))
        nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")
    else:
        ligands = crystal.atoms
        G1 = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table,
                                                     box=crystal.dimensions))
        nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")

    ligand_library = []

    for z, nodes in enumerate(nx.connected_components(G1)):
        if len(ligand_library) == 0:
            ligand_library.append(nodes)
        else:
            Ga = G1.subgraph(nodes)
            is_isomorphic = False
            for ligand_from_library in ligand_library:
                Gb = G1.subgraph(ligand_from_library)

                iso = isomorphism.GraphMatcher(Ga, Gb, node_match=lambda n1, n2: n1['name'] == n2['name'])
                if iso.is_isomorphic():
                    is_isomorphic = True
                    break
            if is_isomorphic == False:
                ligand_library.append(nodes)

    Gtops = [] # nx graphs
    topologies = [] # MDAnalysis universes with the positions
    ligand_tops = [] # parmed topologies

    for idx, nodes in enumerate(ligand_library):
        ligand_coords_mda = MDAnalysis.Merge(crystal.select_atoms(f" index {' '.join(list(map(str, list(nodes)))):}"))

        is_iso, ligand, topology, Gtop = mapping_itp_coords(ligand_coords_mda, ligand_topol) # what if there are more then one topologies (?)
        if not is_iso: # if not isomorphic we paramterize with antechamber:

            ligand_coords_mda.atoms.write('temp.pdb')
            antechamber('temp.pdb', f'linker{idx:}.top')
            ligand = pmd.load_file(f'linker{idx:}.top')

            topology = MDAnalysis.Universe.empty(len(ligand.atoms), trajectory=True)
            topology.add_TopologyAttr('name')
            for a in range(len(topology.atoms)):
                topology.atoms[a].name = ligand.atoms[a].name

            Gtop = nx.Graph([(bond.atom1.idx, bond.atom2.idx) for bond in ligand.bonds])
            nx.set_node_attributes(Gtop, {atom.idx: atom.name[0] for atom in ligand.atoms}, "name")


        ligand_tops.append(ligand)
        topologies.append(topology)
        Gtops.append(Gtop)

    new_ligands = []

    n_ligands = 0

    #parmed has issues if ligands are merged, but molecules are not grouped togheter
    ligand_group = {}

    for z, nodes in enumerate(nx.connected_components(G1)):

        ligands.select_atoms(f" index {' '.join(list(map(str, list(nodes)))):}").write(f"temp.pdb")

        # logger.info(len(nodes))
        Gsub = G1.subgraph(nodes)

        for top_idx, Gtop in enumerate(Gtops):
            new_ligand = topologies[top_idx].copy()

            # Sometimes guesser will not guess right the bonds, decreasing fudge factor might help...
            if len(Gsub) == len(Gtop) and len(Gsub.edges()) != len(Gtop.edges()):
                logger.info("The number of atoms agree... but not the number of bonds, probably bonds are not guessed right")
                fudge_factor = 0.55
                while fudge_factor > 0 and len(Gsub.edges()) != len(Gtop.edges()):
                    Gsub = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(crystal.atoms[Gsub.nodes()],
                                                                             crystal.atoms[Gsub.nodes()].positions,
                                                                             vdwradii=table,
                                                                             box=crystal.dimensions,
                                                                             fudge_factor=fudge_factor))
                    nx.set_node_attributes(Gsub, {atom.index: atom.name[0] for atom in crystal.atoms[Gsub.nodes()]},
                                           "name")
                    fudge_factor -= 0.01
                if fudge_factor < 0.45:

                    raise ValueError("error, cannot guess the bonds")

            iso = isomorphism.GraphMatcher(Gsub, Gtop, node_match=lambda n1, n2: n1['name'] == n2['name'])
            if iso.is_isomorphic():
                for node in list(nodes):
                    # print(node, iso.mapping[node])
                    # print(
                    # print(crystal.atoms[node].name, topology.atoms[iso.mapping[node]].name)
                    new_ligand.atoms[iso.mapping[node]].position = crystal.atoms[node].position
                ligand_group[z] = top_idx

                new_ligands.append(new_ligand.atoms)
                #cage_topol += deepcopy(ligand_tops[top_idx])
                n_ligands+=1
                # print(len(new_ligands),new_ligands )

            else:
                #raise ValueError("Error, cannot match the structures. You need to create topology with nonboned metal manually (or you might try to optimize structure by for example xtb)")
                '''
                print(os.getcwd())
                print(list(nodes))
                #crystal.atoms[list(nodes)].write("temp.pdb")
                #new_ligand.atoms[list(Gsub.nodes())].write("temp.pdb")

                print("The subgraphs are not isomorphic, are you sure these are the same molecules?")
                print(f"Number of atoms: from coordinates: {len(Gsub):d}, from_topology: {len(Gtop):d}")
                print(f"Number of edges: from coordinates: {len(Gsub.edges()):d}, from_topology: {len(Gtop.edges()):d}")

                '''
    cage_topol = None
    if metal_ions is not None:
        for metal_name, indecies in zip(metal_names, metal_indecies):

            metal = create_metal_topol(metal_name, metal_charge, metal_vdw)
            n_metals = len(indecies)
            if cage_topol is None:
                cage_topol = metal * n_metals
            else:
                cage_topol += metal * n_metals


    # we group the ligands togheter
    order = np.argsort(list(ligand_group.values()))

    #firstly the topology file
    for idx in order:
        top_idx = ligand_group[idx]
        cage_topol += deepcopy(ligand_tops[top_idx])

    # and then ligands
    reordered_new_ligands = [metal_ions] + [new_ligands[idx] for idx in order]

    # save new renumered (atoms of linkers and order of linkers) cage
    new_cage = MDAnalysis.Merge(*reordered_new_ligands)
    new_cage.dimensions = crystal.dimensions
    new_cage.atoms.write(output_coord)

    n_metals = len(np.concatenate(metal_indecies))
    cage_topol.write(output_top, [list(range(n_ligands + n_metals))])


    #n_metals = len(metal_indices)

    new_metal_indices = list(range(n_metals))

    return new_metal_indices

import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-metal_name", help="Metal name")
    parser.add_argument("-metal_charge", help="Metal charge")
    parser.add_argument("-output_coord", default="cage.gro", help="Output topology structure of the cage")
    parser.add_argument("-output_top", default="cage.top", help="Output topology structure of the cage")
    parser.add_argument("-metal_vdw", default="merz-tip3p", help="Where to take parameters from")
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    if args.f is not None:
        prepare_initial_topology(args.f, args.metal_name, int(args.metal_charge), args.output_coord, args.output_top, args.metal_vdw)




