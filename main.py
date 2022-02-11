import warnings
warnings.filterwarnings('ignore')

import os
import MDAnalysis
import numpy as np
import networkx as nx

from networkx.algorithms import isomorphism
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array

import argparse
from itertools import permutations
import parmed as pmd
import shutil
from tempfile import mkdtemp

# Taken form: https://gist.github.com/lukasrichters14/
# Dictionary of all elements matched with their atomic masses.
name2mass = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}


class cgbind2gmx():
    path = None
    tmpdir_path = None
    cage = None
    ligand = None
    metal_name = None
    metal_charge =None
    topol_fp = None
    syst_fingerprint = None
    metal_indices =None
    topol_new = None # this is our new topology
    n_metals = None
    n_ligands = None
    name_of_binding_side=None

    output_topol = "cage.top"
    output_coords = "cage.gro"

    def __init__(self, cage_file=None, ligand_file=None, metal_name=None, metal_charge=None, name_of_binding_side=None,
                 cage_cgbind=None, smiles=None, arch_name=None, output_topol=None, output_coords=None):

        print(os.path.dirname(__file__))

        if output_topol is not None:
            self.output_topol = output_topol

        if output_coords is not None:
            self.output_coords = output_coords

        self.tmp_directory()
        print("temp dir", self.tmpdir_path)

    def from_cgbind(self): #TODO
        os.chdir(self.tmpdir_path)

    def from_smiles(self, smiles, arch_name, metal, metal_charge):
        os.chdir(self.tmpdir_path)
        self.create_cage_cgbind(smiles, arch_name, metal, metal_charge)
        self.construct_cage(cage_file='cage.xyz', ligand_file='linker.top', metal_name=metal, metal_charge=metal_charge)
        self.clean_up()

    def from_coords(self, cage_file, ligand_file, metal_name, metal_charge):

        shutil.copy(cage_file, f'{self.tmpdir_path:s}/{cage_file:s}')
        shutil.copy(ligand_file, f'{self.tmpdir_path:s}/{ligand_file:s}')
        os.chdir(self.tmpdir_path)
        self.construct_cage(cage_file=cage_file, ligand_file=ligand_file, metal_name=metal_name, metal_charge=metal_charge)
        self.clean_up()


    def tmp_directory(self):
        self.path = os.getcwd()
        self.tmpdir_path = mkdtemp()


    def clean_up(self):
        # copy everything TODO
        print( f'{self.path:s}/{self.output_coords:s}')
        shutil.copy('cage.gro', f'{self.path:s}/{self.output_coords:s}')
        shutil.copy('cage.top', f'{self.path:s}/{self.output_topol:s}')
        os.chdir(self.path)
        shutil.rmtree(self.tmpdir_path)


    def construct_cage(self, cage_file=None, ligand_file=None, metal_name=None, metal_charge=None, name_of_binding_side=None,):
        self.metal_name = metal_name.upper()
        self.metal_charge = metal_charge
        self.load_cage(cage_file, ligand_file)
        self.load_fingerprint(name_of_binding_side)

        self.prepare_new_topology()

        for metal_index in self.metal_indices:
            mapping_fp_to_new, _ = self.find_mapping(metal_index, self.syst_fingerprint)

            self.adjust_charge(mapping_fp_to_new)

            self.adjust_bonds(mapping_fp_to_new)

            self.adjust_angles(mapping_fp_to_new)

            self.adjust_dihedrals(mapping_fp_to_new)

        self.topol_new.write("cage.top")


    def create_cage_cgbind(self, smiles, arch_name, metal, metal_charge):
        try:
            from cgbind import Linker, Cage
            from antechamber_interface import antechamber
        except:
            raise

        print("[ ] Calling cgbind to create cage")
        linker = Linker(smiles=smiles, arch_name=arch_name)
        cage = Cage(linker, metal=metal, metal_charge=metal_charge)
        cage.print_xyz_file('cage.xyz')

        linker.print_xyz_file('linker.xyz')
        syst = MDAnalysis.Universe('linker.xyz')
        syst.atoms.write('linker.pdb')

        print("[ ] Calling antechamber to parametrize linker") # TODO, that should not be hidden here
        antechamber('linker.pdb', linker.charge, 'linker.top')


    def load_cage(self, cage_file, ligand_file, cutoff=10):
        cage = MDAnalysis.Universe(cage_file)
        cage.atoms.write("temp.gro")
        self.crystal2pdb("temp.gro", ligand_file, "cage.gro", metal_name=self.metal_name)

        self.cage = MDAnalysis.Universe("cage.gro")
        self.ligand = pmd.load_file(ligand_file)

        self.metal_indices = [a for a, name in enumerate(self.cage.atoms.names) if
                         name[:len(self.metal_name)].upper() == self.metal_name]
        self.n_metals = len(self.metal_indices)
        print(self.metal_indices)
        print(self.metal_name)
        print(self.metal_name)
        print("[test]", (len(self.cage.atoms) - self.n_metals) % len(self.ligand.atoms), len(self.cage.atoms),
                                                                                          self.n_metals, len(self.ligand.atoms))
        assert (len(self.cage.atoms) - self.n_metals) % len(self.ligand.atoms) == 0
        self.n_ligands = int((len(self.cage.atoms) - self.n_metals) / len(self.ligand.atoms))

        print(f"Detected M{self.n_metals:d}L{self.n_ligands:d} cage")

        # Find how many ligands are bound to single metal:
        metal_index = self.metal_indices[0] # TODO check if all have the same number of ligands bound to single side

    def find_bound_ligands_nx(self,metal_index, cutoff=10, cutoff_covalent=3):
        cut_sphere = self.cage.select_atoms(f'around {cutoff:f} index {metal_index:d}')
        G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(cut_sphere.atoms, cut_sphere.atoms.positions))
        nx.set_node_attributes(G_cage, {atom.index: atom.name[0] for atom in cut_sphere.atoms}, "name")
        G_sub_cages = [G_cage.subgraph(a) for a in nx.connected_components(G_cage)]

        G_sub_cages_bound = []
        for G_sub_cage in G_sub_cages:
            clusters_of_atoms = self.cage.atoms[G_sub_cage]
            metal = self.cage.atoms[metal_index]
            if np.min(distance_array(clusters_of_atoms.positions, metal.position)) < cutoff_covalent:
                G_sub_cages_bound.append(G_sub_cage)

        return G_sub_cages_bound

    def load_fingerprint(self, name_of_binding_side):
        if self.name_of_binding_side is None:
            self.name_of_binding_side = self.guess_fingerprint()
        #else:
        #    self.name_of_binding_side = name_of_binding_side  # TODO check if it exists

        # Input
        self.topol_fp = pmd.load_file(f'{os.path.dirname(__file__):s}/library/{self.name_of_binding_side:s}.top')
        self.syst_fingerprint = MDAnalysis.Universe(f"{os.path.dirname(__file__):s}/library/{self.name_of_binding_side:s}.pdb")


    def crystal2pdb(self, crystal_pdb, topology_itp, output, metal_name=""):
        crystal = MDAnalysis.Universe(crystal_pdb)
        ligand = pmd.load_file(topology_itp)

        topology = MDAnalysis.Universe.empty(len(ligand.atoms), trajectory=True)
        topology.add_TopologyAttr('name')
        for a in range(len(topology.atoms)):
            topology.atoms[a].name = ligand.atoms[a].name

        table = MDAnalysis.topology.tables.vdwradii
        table['Cl'] = table['CL']

        for atom in crystal.atoms:
            atom.name = atom.name.upper()
        metal_name = metal_name.upper()

        Gtop = nx.Graph([(bond.atom1.idx, bond.atom2.idx) for bond in ligand.bonds])
        nx.set_node_attributes(Gtop, {atom.idx: atom.name[0] for atom in ligand.atoms}, "name")

        if len(metal_name) > 0:
            metal_ions = crystal.select_atoms("name {:s}*".format(metal_name))
            ligands = crystal.select_atoms("not name {:s}*".format(metal_name))
            G1 = nx.Graph(
                MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table))
            nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")
            new_ligands = [metal_ions]
        else:
            ligands = crystal.atoms
            G1 = nx.Graph(
                MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table))
            nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")
            new_ligands = []

        for nodes in nx.connected_components(G1):
            new_ligand = topology.copy()
            # print(len(nodes))
            Gsub = G1.subgraph(nodes)

            iso = isomorphism.GraphMatcher(Gsub, Gtop, node_match=lambda n1, n2: n1['name'] == n2['name'])
            if iso.is_isomorphic():
                for node in list(nodes):
                    # print(node, iso.mapping[node])
                    # print(
                    # print(crystal.atoms[node].name, topology.atoms[iso.mapping[node]].name)
                    new_ligand.atoms[iso.mapping[node]].position = crystal.atoms[node].position

                new_ligands.append(new_ligand.atoms)
                # print(len(new_ligands),new_ligands )
            else:
                print("The subgraphs are not isomorphic, are you sure these are the same molecules?")
                print("Number of atoms:", len(Gsub), len(G1))
                raise Error

        new_cage = MDAnalysis.Merge(*new_ligands)  # , dimensions=crystal.dimensions)
        new_cage.atoms.write(output)

    def guess_fingerprint(self):
        finerprints_names = []
        for file in os.listdir(f"{os.path.dirname(__file__):s}/library"):
            if file.endswith('.pdb'):
                finerprints_names.append(file[:-4])

        metal_index = self.metal_indices[0]

        print("Trying to guess which site it is:", finerprints_names)
        rmsd_best = 1e10
        name_of_binding_side = None
        for finerprint_name in finerprints_names:
            print("[ ] Guessing fingerprint", finerprint_name)
            syst_fingerprint = MDAnalysis.Universe(f"{os.path.dirname(__file__):s}/library/{finerprint_name:s}.pdb")
            mapping_fp_to_new, rmsd = self.find_mapping(metal_index, syst_fingerprint, guessing=True)
            if rmsd < rmsd_best:
                rmsd_best = rmsd
                name_of_binding_side = finerprint_name

        print("[+] Best fingerprint", name_of_binding_side, "rmsd", rmsd_best)
        if rmsd_best > 1.0:
            print("[!] Rmsd is quite large, want to proceed?") #TODO
        return name_of_binding_side

    #TODO assert rediculsy small ligands

    def find_mapping(self, metal_index, syst_fingerprint, cutoff=9, guessing=False):
        G_fingerprint = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(syst_fingerprint.atoms, syst_fingerprint.atoms.positions))
        nx.set_node_attributes(G_fingerprint, {atom.index: atom.name[0] for atom in syst_fingerprint.atoms}, "name")
        G_fingerprint_subs = [G_fingerprint.subgraph(a) for a in nx.connected_components(G_fingerprint)]

        print("     [ ] Mapping fingerprint to metal center:", metal_index)
        '''
        cut_sphere = self.cage.select_atoms(f'index {metal_index:d} or around {cutoff:f} index {metal_index:d}')
        G_cage = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(cut_sphere.atoms, cut_sphere.atoms.positions))
        nx.set_node_attributes(G_cage, {atom.index: atom.name[0] for atom in cut_sphere.atoms}, "name")
        G_sub_cages = [G_cage.subgraph(a) for a in nx.connected_components(G_cage)]
        G_sub_cages = sorted(G_sub_cages, key=len, reverse=True) # we assume that the largerst group are  ligands, is this reasonable? TODO
        '''

        G_sub_cages = self.find_bound_ligands_nx(metal_index)
        number_ligands_bound = len(G_sub_cages)

        print("         Number of ligands bound to metal", number_ligands_bound)

        if len(G_sub_cages) != len(G_fingerprint_subs) and guessing:
            print("[!] Not the same number of sites", guessing)
            return None, 1e10
        elif len(G_sub_cages) != len(G_fingerprint_subs) and not guessing:
            print("[!] Not the same number of sites", guessing)
            raise

        selected_atoms = [metal_index]
        for G_sub_cage in G_sub_cages:

            if guessing:  # if we want to guess site, that might be not fulfiled, and that means it is not the corret site
                if not len(G_fingerprint_subs[0]) < len(G_sub_cage):
                    print("not subgroup", len(G_fingerprint_subs[0]), len(G_sub_cage))
                    return None, 1e10
            else:
                # Check if cutted ligands from cage are larger then the ligands in finerprint (they should be if cutoff is 10!)
                assert len(G_fingerprint_subs[0]) < len(G_sub_cage)

            ismags = nx.isomorphism.ISMAGS(G_sub_cage, G_fingerprint_subs[0],
                                           node_match=lambda n1, n2: n1['name'] == n2['name'])

            largest_common_subgraph = list(ismags.largest_common_subgraph())

            if guessing:  # if we want to guess site, that might be not fulfiled, and that means it is not the corret site
                if not len(G_fingerprint_subs[0]) == len(largest_common_subgraph[0]):
                    return None, 1e10
            else:
                # Check if the subgraph is the same size as fingerprint,
                assert len(G_fingerprint_subs[0]) == len(largest_common_subgraph[0])

            selected_atoms += largest_common_subgraph[0].keys()

        cut_sphere = self.cage.select_atoms(f'index {metal_index:d} or around {cutoff:f} index {metal_index:d}')
        connected_cut_system = cut_sphere.select_atoms("index " + " ".join(map(str, selected_atoms)))

        G_site = nx.Graph(
            MDAnalysis.topology.guessers.guess_bonds(connected_cut_system.atoms, connected_cut_system.atoms.positions))
        nx.set_node_attributes(G_site, {atom.index: atom.name[0] for atom in connected_cut_system.atoms}, "name")
        G_site_subs = [G_site.subgraph(a) for a in nx.connected_components(G_site)]

        best_rmsd = 1e10
        best_mapping = None

        indecies2number = {}
        for a, b in enumerate(connected_cut_system.indices):
            indecies2number[b] = a

        permutations_list = permutations(range(number_ligands_bound))

        for perm in list(permutations_list):
            mappings = []
            for a, b in enumerate(perm):
                iso = isomorphism.GraphMatcher(G_fingerprint_subs[a], G_site_subs[b],
                                               node_match=lambda n1, n2: n1['name'] == n2['name'])
                mappings.append([subgraph_mapping for subgraph_mapping in iso.subgraph_isomorphisms_iter()])

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

                # TODO REMOVE METAL,or ADD ITS INDEX!
                reordered_system = connected_cut_system.atoms[
                    [0] + [indecies2number[value] for value in new_new_dict.values()]]  # Find where metal is
                _, rmsd = align.alignto(syst_fingerprint, reordered_system)

                # print(rmsd)
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_mapping = new_new_dict

        temp = np.where([syst_fingerprint.atoms.names == self.metal_name])[0]
        if len(temp) == 1:
            best_mapping[temp[0]] = metal_index
        else:
            print("ERROR, more than one metal in the site")
            raise

        return best_mapping, best_rmsd

    def prepare_new_topology(self):
        metal = pmd.load_file(f'{os.path.dirname(__file__):s}/library/M.itp')
        metal.atoms[0].name = self.metal_name
        metal.atoms[0].charge = self.metal_charge
        metal.atoms[0].mass = name2mass[self.metal_name]

        cage_topol = metal * self.n_metals + self.ligand * self.n_ligands
        cage_topol.write('test_cage.top', [list(range(self.n_ligands + self.n_metals))])

        self.topol_new = pmd.load_file('test_cage.top', parametrize=False)  # , skip_bonds=True)
        # we need to copy paramters back, for some reasons if paramtrize is False
        # parmed is sometimes wierd :-(

        topol_new2 = pmd.load_file('test_cage.top')
        for a in range(len(self.topol_new.atoms)):
            self.topol_new.atoms[a].type = topol_new2.atoms[a].type
            self.topol_new.atoms[a].epsilon = topol_new2.atoms[a].epsilon
            self.topol_new.atoms[a].sigma = topol_new2.atoms[a].sigma
            self.topol_new.atoms[a].rmin = topol_new2.atoms[a].rmin
            self.topol_new.atoms[a].charge = topol_new2.atoms[a].charge

    def adjust_charge(self, mapping_fp_to_new):
        print("   [ ] Changing charges and atomtypes")
        sum_of_charge_diffrences = 0

        for a in mapping_fp_to_new:
            print("         ", self.topol_new.atoms[mapping_fp_to_new[a]].type, self.topol_new.atoms[mapping_fp_to_new[a]].name,
                  self.topol_new.atoms[mapping_fp_to_new[a]].charge, "-->",
                  self.topol_fp.atoms[a].type, self.topol_fp.atoms[a].name, self.topol_fp.atoms[a].charge)
            atom = self.topol_fp.atoms[a]

            if atom.type not in self.topol_new.parameterset.atom_types.keys():
                print("      [^] Adding new atomtype:", atom.type)
                atomtype = pmd.topologyobjects.AtomType(name=atom.type, number=atom.number, mass=atom.mass)
                self.topol_new.parameterset.atom_types[atom.type] = atomtype

            self.topol_new.atoms[mapping_fp_to_new[a]].type = self.topol_fp.atoms[a].type
            self.topol_new.atoms[mapping_fp_to_new[a]].epsilon = self.topol_fp.atoms[a].epsilon
            self.topol_new.atoms[mapping_fp_to_new[a]].sigma = self.topol_fp.atoms[a].sigma
            self.topol_new.atoms[mapping_fp_to_new[a]].rmin = self.topol_fp.atoms[a].rmin
            self.topol_new.atoms[mapping_fp_to_new[a]].charge += self.topol_fp.atoms[a].charge  # topol_fp.atoms[a].charge

            sum_of_charge_diffrences += self.topol_fp.atoms[a].charge

    def adjust_bonds(self, mapping_fp_to_new):
        print("   [ ] Adding new bonds to topology")
        for bond_fp in self.topol_fp.bonds:
            found = False
            for bond_new in self.topol_new.bonds:

                # Check if atoms in mapping (that they are not the additional atoms)
                if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                    if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                         mapping_fp_to_new[bond_fp.atom2.idx]) or
                            ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                              mapping_fp_to_new[bond_fp.atom1.idx]))):
                        if (bond_fp.type != bond_new.type):
                            print("      [o] Diffrent bond type ", bond_new.type, bond_fp.type)
                            type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                                          list=self.topol_new.bond_types)
                            # if type_to_assign not in topol_new.bond_types:
                            self.topol_new.bond_types.append(type_to_assign)
                            bond_new.type = bond_fp.type

        print("   [ ] Adding new bonds to topology")

        for bond_fp in self.topol_fp.bonds:
            found = False

            for bond_new in self.topol_new.bonds:
                # if bond_fp.atom1.name=="ZN":
                #    print(bond_fp.atom1.idx , bond_new)

                # print(bond_new, bond_new.atom1.name)
                # Check if atoms in mapping (that they are not the additional atoms)
                if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                    # print("It is in the mapping")

                    if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                         mapping_fp_to_new[bond_fp.atom2.idx]) or
                            ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                              mapping_fp_to_new[bond_fp.atom1.idx]))):
                        if (bond_fp.type != bond_new.type):
                            print("      [o] Diffrent bond type ", bond_new.type, bond_fp.type)
                        found = True
                else:
                    found = True

            if not found:
                # type_to_assign= bond_fp.type
                # type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq)
                # print("AAA", bond_fp.type.k, bond_fp.type.req, len(topol_new.bond_types))

                type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                              list=self.topol_new.bond_types)

                print("      [o] New bond:", mapping_fp_to_new[bond_fp.atom1.idx],
                      mapping_fp_to_new[bond_fp.atom2.idx], type_to_assign)
                # if type_to_assign not in topol_new.bond_types:
                self.topol_new.bond_types.append(type_to_assign)

                atom1 = self.topol_new.atoms[mapping_fp_to_new[bond_fp.atom1.idx]]
                atom2 = self.topol_new.atoms[mapping_fp_to_new[bond_fp.atom2.idx]]
                self.topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))


    def adjust_angles(self, mapping_fp_to_new):
        print("   [ ] Adding new angles to topology")
        for angle_fp in self.topol_fp.angles:
            found = False
            if angle_fp.atom1.idx in mapping_fp_to_new and angle_fp.atom2.idx in mapping_fp_to_new and angle_fp.atom3.idx in mapping_fp_to_new:
                for angle_new in self.topol_new.angles:
                    #Check if atoms in mapping (that they are not the additional atoms)
                    if ((angle_new.atom1.idx==mapping_fp_to_new[angle_fp.atom1.idx] and angle_new.atom2.idx==mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx==mapping_fp_to_new[angle_fp.atom3.idx]) or
                        ((angle_new.atom1.idx==mapping_fp_to_new[angle_fp.atom3.idx] and angle_new.atom2.idx==mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx==mapping_fp_to_new[angle_fp.atom1.idx]))):
                        if (angle_fp.type!=angle_new.type):
                            print("      [o] Diffrent angle type ",angle_new.type,"->", angle_fp.type )
                            #angle_new.funct = angle_fp.funct
                            #angle_new.type.k = angle_fp.type.k
                            #angle_new.type.theteq = angle_fp.type.theteq


                            type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq, list = self.topol_new.angle_types)

                            #if type_to_assign not in topol_new.angle_types:
                            self.topol_new.angle_types.append(type_to_assign)
                            angle_new.type = type_to_assign

                        found = True
            else:
                    print("[-] Not found:", angle_fp.atom1.idx+1, angle_fp.atom2.idx+1, angle_fp.atom3.idx+1 )
                    print("              ", self.topol_fp.atoms[angle_fp.atom1.idx])
                    print("              ", self.topol_fp.atoms[angle_fp.atom2.idx])
                    print("              ", self.topol_fp.atoms[angle_fp.atom3.idx])
                    print("And you should be worry if missing")
                    found = True
                    raise

            if not found:
                #type_to_assign= angle_fp.type
                type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq,  list = self.topol_new.angle_types )
                #if type_to_assign not in topol_new.angle_types:
                self.topol_new.angle_types.append(type_to_assign)
                print("new type", len(self.topol_new.angle_types))

                print("      [x] New angle:",angle_fp.atom1.name,'-',angle_fp.atom2.name,'-',angle_fp.atom3.name, mapping_fp_to_new[angle_fp.atom1.idx], mapping_fp_to_new[angle_fp.atom2.idx], mapping_fp_to_new[angle_fp.atom3.idx], type_to_assign)
                #print("              ", topol_fp.atoms[angle_fp.atom1.idx])
                #print("              ", topol_fp.atoms[angle_fp.atom2.idx])
                #print("              ", topol_fp.atoms[angle_fp.atom3.idx])
                atom1= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom1.idx]]
                atom2= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom2.idx]]
                atom3= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom3.idx]]

                self.topol_new.angles.append(pmd.topologyobjects.Angle(atom1, atom2, atom3, type=type_to_assign))

    def adjust_dihedrals(self, mapping_fp_to_new):
        print("   [ ] Adding new dihedrals to topology")
        for dihedral_fp in self.topol_fp.dihedrals:
            found = False
            if dihedral_fp.atom1.idx in mapping_fp_to_new and dihedral_fp.atom2.idx in mapping_fp_to_new and dihedral_fp.atom3.idx in mapping_fp_to_new and dihedral_fp.atom4.idx in mapping_fp_to_new:
                for dihedral_new in self.topol_new.dihedrals:
                    #Check if atoms in mapping (that they are not the additional atoms)
                    if (((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom1.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom4.idx]) or
                        ((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom4.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom1.idx]))) and
                       (dihedral_new.type.per==dihedral_fp.type.per)):
                        if (dihedral_fp.type!=dihedral_new.type):
                            print("      [a] Diffrent Dihedral type ",dihedral_fp.atom1.name,"-",dihedral_fp.atom2.name,"-",dihedral_fp.atom3.name,"-",dihedral_fp.atom4.name)
                            print("          ", dihedral_new.type)
                            print("        ->", dihedral_fp.type)

                            type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedral_fp.type.phi_k, per=dihedral_fp.type.per,
                                                                             phase=dihedral_fp.type.phase, scee=dihedral_fp.type.scee,
                                                                             scnb=dihedral_fp.type.scnb, list = self.topol_new.dihedral_types)
                            #if type_to_assign not in topol_new.dihedral_types:
                            self.topol_new.dihedral_types.append(type_to_assign)
                            dihedral_new.type = type_to_assign

                        found = True
            else:
                    print("[-] Not found:", dihedral_fp.atom1.idx+1, dihedral_fp.atom2.idx+1, dihedral_fp.atom3.idx+1 , dihedral_fp.atom4.idx+1 )
                    print("              ", self.topol_fp.atoms[dihedral_fp.atom1.idx])
                    print("              ", self.topol_fp.atoms[dihedral_fp.atom2.idx])
                    print("              ", self.topol_fp.atoms[dihedral_fp.atom3.idx])
                    print("              ", self.topol_fp.atoms[dihedral_fp.atom4.idx])
                    print("And you should be worry if missing")
                    found = True
                    raise

            if not found:
                type_to_assign= pmd.topologyobjects.DihedralType(phi_k=dihedral_fp.type.phi_k, per=dihedral_fp.type.per,
                                                                phase=dihedral_fp.type.phase, scee=dihedral_fp.type.scee,
                                                                scnb=dihedral_fp.type.scnb, list = self.topol_new.dihedral_types)
                #if type_to_assign not in topol_new.dihedral_types:
                self.topol_new.dihedral_types.append(type_to_assign)
                print("      [b] New dihedral:",dihedral_fp.atom1.name,"-",dihedral_fp.atom2.name, mapping_fp_to_new[dihedral_fp.atom1.idx], mapping_fp_to_new[dihedral_fp.atom2.idx], mapping_fp_to_new[dihedral_fp.atom3.idx], type_to_assign)

                atom1 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom1.idx]]
                atom2 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom2.idx]]
                atom3 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom3.idx]]
                atom4 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom4.idx]]

                self.topol_new.dihedrals.append(pmd.topologyobjects.Dihedral(atom1, atom2, atom3, atom4, type=type_to_assign))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Structure pf the cage ")
    parser.add_argument("-linker_topol", help="Topology of the linker ") # TODO
    parser.add_argument("-linker_coords", help="Structure of the linker ")  # TODO
    parser.add_argument("-fingerprint", default=None, help="Structure of the linker ")  # TODO

    parser.add_argument("-name", default='UNK', help="trajectory (traj_comp.xtc)")
    parser.add_argument("-smiles", default=None, help="Smiles of the linker")
    parser.add_argument("-arch_name", default=None, help="Architecture")
    parser.add_argument("-metal", default='Pd', help="Metal name")
    parser.add_argument("-metal_charge", default='2', help="Metal charge")

    parser.add_argument("-o", default='cage.gro', help="Output coordination file")
    parser.add_argument("-ot", default='cage.top', help="Output topology file")

    return parser.parse_args()

# MAKE it less louder
if __name__ == '__main__':
    args = get_args()
    cgbind2gmx = cgbind2gmx(output_coords=args.o, output_topol=args.ot)

    if args.fingerprint is not None:
        cgbind2gmx.name_of_binding_side = args.fingerprint

    if args.smiles is not None and args.arch_name is not None and args.metal is not None and args.metal_charge is not None:
        cgbind2gmx.from_smiles(args.smiles, args.arch_name, args.metal, int(args.metal_charge))
    elif args.f is not None and args.linker_topol is not None and args.metal is not None and args.metal_charge is not None:
        cgbind2gmx.from_coords(args.f, args.linker_topol, args.metal, int(args.metal_charge))
    else:
        print("One of the three values needs to be specified:")
        print('a) smiles, arch_name, metal, charge')
        print("b) f, linker_topol, metal, charge")


    #cgbind2gmx(cage_file="/media/piskorz/Oxford/Zincage/big_cage/cage.pdb", ligand_file='/media/piskorz/Oxford/Zincage/parameters/large_ligand/ZnL.top',
    #           metal_name = "ZN", metal_charge = 2, name_of_binding_side = "ZnN2")
