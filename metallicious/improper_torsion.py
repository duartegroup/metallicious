import MDAnalysis
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np

import os
from scipy.optimize import minimize
import networkx as nx
from metallicious.mapping import map_two_structures
from metallicious.utils import new_directory, strip_numbers_from_atom_names
from metallicious.log import logger

def find_donor_indices(bonds):
    """
    Returns indices of atoms connected to metal (assumes metal has index 0)
    :param bonds: (list(int,int)) list of indices pairs describing bonds
    :return: (list(int)) indices of donor atoms
    """
    indices = list(nx.generators.ego_graph(nx.Graph([bond for bond in bonds]), 0, radius=1).nodes)
    indices.remove(0)

    return indices


def find_potential_impropers(bonds, donor_indices):
    """
    Finds improper dihedral: it assumes that donor atom is 'out-of-plane' and therefore 3 atoms which it connects
    (including metal) should form improper dihedral
    :param bonds: (list(int,int)) list of indices pairs
    :param donor_indices: (list(int)) list of indices
    :return: list(list(int,int,int,int)) list of dihedrals
    """
    indices = []
    for donor_index in donor_indices:
        indices_connected_to_donor = list(
            nx.generators.ego_graph(nx.Graph([bond for bond in bonds]), donor_index, radius=1).nodes)
        if len(indices_connected_to_donor) == 4:
            indices_of_middle_atoms = [index for index in indices_connected_to_donor if index not in [0, donor_index]]
            indices.append([0, indices_of_middle_atoms[0], indices_of_middle_atoms[1], donor_index])
    return indices

def filter_aromatic_impropers(filename, potential_impropers, metal_name):
    """
    Filters found potential dihedrals for the on which consist of aromatic subunits
    :param filename: (str) filename of coordination file
    :param potential_impropers: (list(list(int,int,int,int)) list of potential dihedrals
    :param metal_name: (str) name of file
    :return:(list(list(int,int,int,int)) list of aromatic dihedrals
    """
    aromatic_impropers = []

    syst = MDAnalysis.Universe(filename)
    for atom in syst.select_atoms(f"name {metal_name:}"):
        atom.type = 'H'  # we change type of metal for hydrogen, because VdW radii is required for aromaticity

    if not hasattr(syst.atoms[0], 'element'):
        guessed_elements = MDAnalysis.topology.guessers.guess_types(strip_numbers_from_atom_names(syst.names))
        for atom in syst.select_atoms(f"name {metal_name:}"):
            guessed_elements[atom.idx] = 'H'
        syst.universe.add_TopologyAttr('elements', guessed_elements)

    aromaticity = MDAnalysis.topology.guessers.guess_aromaticities(syst.atoms)

    for improper in potential_impropers:
        if aromaticity[improper[1:]].all():
            aromatic_impropers.append(improper)
    return aromatic_impropers


def scan_improper(improper, filename='bonded/site0_optimised_orca.xyz', charge=0, mult=1, x=np.linspace(0, 5, 2)):
    """
    Modification to autodE, which allows to perform scan of the dihedrals
    :param improper: list(int,int,int,int) dihedral for QM scan
    :param filename: (str) filename of the input structure
    :param charge: (int) charge for QM calculations
    :param mult: (int) multiplicity for QM calculations
    :param x: list(float): space which will be scanned
    :return: list(float): electronic energy of the scanned coordinate
    """
    import autode as ade
    from collections.abc import MutableMapping
    from autode.wrappers.ORCA import logger, print_added_internals, print_distance_constraints, \
        print_cartesian_constraints, print_num_optimisation_steps, print_point_charges, print_default_params, \
        print_coordinates
    from autode.values import Angle

    class DihedralConstraints(MutableMapping):
        def __init__(self, *args, **kwargs):
            self._store = dict()
            self.update(dict(*args, **kwargs))  # use the free update to set keys

        def __getitem__(self, key):
            return self._store[self._key_transform(key)]

        def __delitem__(self, key):
            del self._store[self._key_transform(key)]

        def __iter__(self):
            return iter(self._store)

        def __len__(self):
            return len(self._store)

        @staticmethod
        def _key_transform(key):
            """Transform the key to a sorted tuple"""
            return tuple(key)

        def __setitem__(self, key, value):
            """
            Set a key-value pair in the dictionary
            -----------------------------------------------------------------------
            Arguments:
                key (tuple(int)): Pair of atom indexes
                value (int | float): Distance
            """
            try:
                n_unique_atoms = len(set(key))
            except TypeError:
                raise ValueError(f"Cannot set a key with {key}, must be iterable")
            if n_unique_atoms != 4:
                logger.error(
                    "Tried to set a distance constraint with a key: "
                    f"{key}. Must be a unique pair of atom indexes"
                )
                return
            if any(int(atom_idx) < 0 for atom_idx in key):
                raise ValueError(
                    "Distance constraint key must be an atom index "
                    f"pair but had: {key} which cannot be valid (<0)"
                )
            self._store[self._key_transform(key)] = Angle(value)

    def print_dihedral_constraints(inp_file, molecule):
        """Print the distance constraints to the input file"""
        logger.debug("Adding dihedral to file")
        if molecule.constraints.dihedral is None:
            return
        print("%geom Constraints", file=inp_file)
        for (i, j, k, l), ang in molecule.constraints.dihedral.items():
            print("{ D", i, j, k, l, ang, "C }", file=inp_file)
        print("    end\nend", file=inp_file)
        return

    def generate_input_for_new(self, calc: "CalculationExecutor") -> None:
        logger.debug("Using new input file")
        molecule = calc.molecule
        keywords = self.get_keywords(calc.input, molecule)
        # with open(calc.input.directory, "w") as inp_file: #change from 1.2.1 autode to 1.3
        with open(calc.input.filename, "w") as inp_file:
            print("!", *keywords, file=inp_file)
            self.print_solvent(inp_file, molecule, keywords)
            print_added_internals(inp_file, calc.input)
            print_distance_constraints(inp_file, molecule)
            print_dihedral_constraints(inp_file, molecule)
            print_cartesian_constraints(inp_file, molecule)
            print_num_optimisation_steps(inp_file, molecule, calc.input)
            print_point_charges(inp_file, calc.input)
            print_default_params(inp_file)
            if calc.n_cores > 1:
                print(f"%pal nprocs {calc.n_cores}\nend", file=inp_file)
            print_coordinates(inp_file, molecule)
        return None

    def generate_input_for(self, calc: "CalculationExecutor") -> None:
        molecule = calc.molecule
        keywords = self.get_keywords(calc.input, molecule)

        # with open(calc.input.directory, "w") as inp_file:
        with open(calc.input.filename, "w") as inp_file:
            print("!", *keywords, file=inp_file)
            self.print_solvent(inp_file, molecule, keywords)
            print_added_internals(inp_file, calc.input)
            print_distance_constraints(inp_file, molecule)
            print_cartesian_constraints(inp_file, molecule)
            print_num_optimisation_steps(inp_file, molecule, calc.input)
            print_point_charges(inp_file, calc.input)
            print_default_params(inp_file)

            if calc.n_cores > 1:
                print(f"%pal nprocs {calc.n_cores}\nend", file=inp_file)

            print_coordinates(inp_file, molecule)
        return None

    method = ade.methods.ORCA()
    if method.is_available is False:
        raise NameError("For parametrization of templates, QM software ORCA is required")

    # overwrite the autode ORCA input
    ade.wrappers.ORCA.ORCA.generate_input_for = generate_input_for_new

    molecule = ade.Molecule(filename, charge=charge, mult=mult)
    name = molecule.name
    energies = []

    here = os.getcwd()

    dir_name = f"improper_{'_'.join(map(str, improper)):s}"
    new_directory(dir_name)
    os.chdir(dir_name)

    for angle in x:
        logger.debug(f"{molecule.atoms[0]:}, angle: {angle}")
        molecule.name = name + '_' + str(angle)
        molecule.constraints.dihedral = DihedralConstraints({tuple(improper): angle})
        molecule.optimise(method=method)
        energies.append(molecule.energy.to("kcal"))  # amber uses kcal/mol units
    energies = np.array(energies) - np.min(energies)  # we change it to array and normalize it

    logger.info(f"Torsion energy: {energies}")

    # we reverse autodE orca input
    ade.wrappers.ORCA.ORCA.generate_input_for = generate_input_for
    os.chdir(here)
    return energies


def evaluate_angle(improper, filename='bonded/site0_optimised_orca.xyz'):
    """
    Calculates dihedral angle for input structure
    :param improper: list(int) list of indices indicating dihedral
    :param filename: (str) filename
    :return: (float) dihedral angle (in degrees)
    """

    syst = MDAnalysis.Universe(filename)
    positions = syst.atoms[improper].positions
    angle = np.rad2deg(calc_dihedrals(*positions))
    return angle


def evaluate_improper_energy(energies, x=np.linspace(0, 20, 5)):
    """
    Calculate force constant assuming harmonic force
    :param energies: (list(float)) input energies
    :param x: (list(float)) values of scanned coordinate
    :return: (float) force-constant
    """

    initial_force = 0

    def evaluate(k):
        return np.sum(np.abs(k * (np.deg2rad(x)) ** 2 - energies))  # in AMBER format, that is why it's missing 0.5

    result = minimize(evaluate, x0=[initial_force], options={'eps': 1})
    return result.x[0]


def improper_value_calculation(improper, filename='bonded/site_optimised_orca.xyz', charge=0, mult=1,
                               angles_to_check=np.linspace(0, 20, 5), angle_cutoff=5):
    """
    For given improper dihedral perform QM scan

    :param improper: (list(int)) indices of atoms forming improper dihedral
    :param filename: (str) filename for the input structure
    :param charge: (int) charge for QM calculations
    :param mult: (int) multiplicity for QM calculations
    :param angles_to_check: (list(float)) space to scan
    :param angle_cutoff: (float) consider only dihedral which form angle smaller than angle_cutoff
    :return:
    """
    angle = evaluate_angle(improper, filename=filename)
    if angle < angle_cutoff:
        angle = 0

        energies = scan_improper(improper, filename=filename, charge=charge, mult=mult, x=angles_to_check)
        value = evaluate_improper_energy(energies, x=angles_to_check)

        logger.info("\tEnergies:", energies)
        logger.info("\tForce constant", value)

        return angle, value
    else:
        return angle, 0.0


def find_impropers_connected_to_metal(bonds, metal_name, filename='bonded/site_opt_orca.xyz'):
    """
    Finds improper dihedrals connected to metal and aromatic moiety
    :param bonds: (list(int,int)) list of indices pairs describing bonds
    :param metal_name: (str) metal name
    :param filename: (str) filename of the input structure
    :return: (list(list(int)) list of dihedrals
    """
    donor_indices = find_donor_indices(bonds)
    potential_impropers = find_potential_impropers(bonds, donor_indices)
    impropers = filter_aromatic_impropers(filename, potential_impropers, metal_name)
    return impropers


def symmetric_improper(improper):
    """
    For given improper dihedral give equivalent dihedral but with reordered atoms
    :param improper: list(int) indices of atoms forming dihedral
    :return: list(int) indices of atoms forming dihedral
    """
    return [improper[0], improper[2], improper[1], improper[3]]


def find_impropers_and_values(bonds, metal_name, unique_ligands_pattern, starting_index, indices, charge, mult=1,
                              filename='bonded/site_opt_orca.xyz'):
    """
    The main routine for identifying dihedrals

    :param bonds: (list(int)) list of indices pairs describing bonds
    :param metal_name: (str) name of metal
    :param unique_ligands_pattern: (list(int)) list describing which ligands are the same
    :param starting_index: (list(int)) indices of the first atoms in residues
    :param indices: (list(list(int)) list of ligand indices
    :param charge: (int) charge of the system
    :param mult: (int) multiplicity of the system
    :param filename: (str) filename of the system
    :return: dict() describing dihedral {(id1, id2, id3, id4):(angle,force_constant)}
    """

    logger.info("Calculating improper dihedrals")
    impropers = find_impropers_connected_to_metal(bonds, metal_name, filename=filename)

    for unique_ligand in list(set(unique_ligands_pattern)):
        new_dihedrals = {}

        impropers_to_calculate = []

        same_ligands = np.where(unique_ligand == np.array(unique_ligands_pattern))[0]

        ligand_indices = indices[1:][same_ligands[0]]
        first_start_index = starting_index[1]

        for improper in impropers:
            if sum([1 for a in ligand_indices if a in improper]) == 3:
                impropers_to_calculate.append(improper)
        logger.debug(f"Improper to calculate: {impropers_to_calculate}")

        # calculating
        values_of_checked_impropers = []
        for improper in impropers_to_calculate:
            value = improper_value_calculation(improper, filename=filename, charge=charge, mult=mult)

            values_of_checked_impropers.append(value)
            new_dihedrals[tuple(improper)] = tuple(value)
            # new_dihedrals[tuple([improper[0],improper[2],improper[1],improper[3]])] = tuple(value)

        if len(same_ligands) > 1:
            for same_ligand in same_ligands[1:]:
                logger.debug(f"{same_ligand:}, {starting_index[1:][same_ligand] - first_start_index:}")

                syst = MDAnalysis.Universe(filename)
                ligand1 = syst.atoms[
                    [0] + indices[1:][same_ligands[0]]]  # ligand + metal (metal is necessary for symmetry)
                ligand2 = syst.atoms[[0] + indices[1:][same_ligand]]

                mapping, _ = map_two_structures(0, ligand2, ligand1,
                                                metal_name)  # we map one of the ligands+metal on second
                proposed_torsions = [[mapping[a] for a in improper] for improper in
                                     impropers_to_calculate]  # map calculated impropers on the mapped ligands

                # check if found the same torsions:
                assert sum(1 for proposed_torsion in proposed_torsions if (list(proposed_torsion) in impropers) or (
                    symmetric_improper(proposed_torsion)) in impropers) == len(proposed_torsions)
                for idx_tor, proposed_torsion in enumerate(proposed_torsions):
                    new_dihedrals[tuple(proposed_torsion)] = tuple(values_of_checked_impropers[idx_tor])

    return new_dihedrals
