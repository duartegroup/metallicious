
import MDAnalysis
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from collections.abc import MutableMapping
from autode.values import Angle
import autode as ade
from scipy.optimize import minimize

method = ade.methods.ORCA()

def find_donor_indices(bonds_with_names, metal_name):
    indicies = []
    for bond_name in bonds_with_names:
        if metal_name.title() in bond_name[1]:
            if bond_name[1][0] == metal_name.title():
                indicies.append(bond_name[0][1])
            elif bond_name[1][1] == metal_name.title():
                indicies.append(bond_name[0][0])
    return indicies


def find_potential_impropers(bonds_with_names, donor_indices, metal_name):
    indices = []
    for donor_index in donor_indices:

        indices.append([donor_index])
        metal_index = None
        for bond_name in bonds_with_names:
            if donor_index in bond_name[0]:
                if donor_index == bond_name[0][0]:
                    if bond_name[1][1] == metal_name.title():  # we check if not metal, we want to add metal at the end
                        metal_index = bond_name[0][1]
                    else:
                        indices[-1].append(bond_name[0][1])

                elif donor_index == bond_name[0][1]:
                    if bond_name[1][0] == metal_name.title():  # we check if not metal, we want to add metal at the end
                        metal_index = bond_name[0][0]
                    else:
                        indices[-1].append(bond_name[0][0])
        indices[-1].append(metal_index)
    return indices


def filter_aromatic_impropers(filename, potential_impropers, metal_name):
    print(filename, potential_impropers, metal_name)
    aromatic_impropers = []

    syst = MDAnalysis.Universe(filename)
    for atom in syst.select_atoms(f"name {metal_name:}"):
        atom.type = 'H'  # we change type of metal for hydrogen, becasue vdwradii is required for aromaticity

    if not hasattr(syst.atoms[0], 'element'):
        guessed_elements = MDAnalysis.topology.guessers.guess_types(syst.names)
        for atom in syst.select_atoms(f"name {metal_name:}"):
            guessed_elements[atom.idx] = 'H'
        syst.universe.add_TopologyAttr('elements', guessed_elements)

    aromaticity = MDAnalysis.topology.guessers.guess_aromaticities(syst.atoms)
    # print(aromaticity)
    for improper in potential_impropers:
        if aromaticity[improper[:-1]].all():
            aromatic_impropers.append(improper)
    return aromatic_impropers


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
    if molecule.constraints.dihedral is None:
        return
    print("%geom Constraints", file=inp_file)
    for (i, j, k, l), ang in molecule.constraints.dihedral.items():
        print("{ D", i, j, k, l, ang, "C }", file=inp_file)
    print("    end\nend", file=inp_file)
    return


from autode.wrappers.ORCA import *

def generate_input_for_new(self, calc: "CalculationExecutor") -> None:
    molecule = calc.molecule
    keywords = self.get_keywords(calc.input, molecule)
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

def scan_improper(improper, filename='bonded/site0_optimised_orca.xyz', charge=0, x=np.linspace(0,30, 7)):
    #overwrite the autode ORCA input
    ade.wrappers.ORCA.ORCA.generate_input_for = generate_input_for_new

    molecule = ade.Molecule(filename, charge = charge)
    name = molecule.name
    energies = []
    for angle in x:
        print(molecule.atoms[0], "angle", angle)
        molecule.name = name + '_' + str(angle)
        molecule.constraints.dihedral = DihedralConstraints({tuple(improper):angle})
        molecule.optimise(method=method)
        energies.append(molecule.energy.to("kcal")) #amber uses kcal/mol units
    energies = np.array(energies) - np.min(energies) # we change it to array and normalize it

    print("Trosion energy:", energies)
    # we reverse autode orca input
    ade.wrappers.ORCA.ORCA.generate_input_for = generate_input_for

    return energies


def evaluate_angle(improper, filename='bonded/site0_optimised_orca.xyz'):
    syst = MDAnalysis.Universe(filename)
    positions = syst.atoms[improper].positions
    angle = np.rad2deg(calc_dihedrals(*positions))
    if np.abs(angle)<10:
        return 0
    else:
        print(f"Unusall value of improper angle: {angle:}")
        raise

def evalulate_improper_energy(energies, x=np.linspace(0,30, 7)):
    initial_force = 0
    def evaluate(k):
        return np.sum(np.abs(k*(np.deg2rad(x))**2-energies)) # This is in AMBER format, that is why it's missing 0.5
    result = minimize(evaluate, x0=[initial_force], options = {'eps': 1})
    return result.x[0]


def improper_value_calculation(improper, filename='bonded/site0_optimised_orca.xyz', charge=0):
    angle = evaluate_angle(improper, filename=filename)
    energies = scan_improper(improper, filename=filename, charge=charge, x=np.linspace(0,30, 7))
    value = evalulate_improper_energy(energies, x=np.linspace(0,30, 7))
    print("Force constant", value)

    value *=0.5 # we devide by half, becasue we will double count it, with two improper dihedrals for symmetry
    return (angle, value)

def dummy_calculation(improper):
    return (0,10)

def find_impropers_connected_to_metal(bonds_with_names, metal_name, filename='bonded/site0_opt_orca.xyz'):
    donor_indices = find_donor_indices(bonds_with_names, metal_name)
    potential_impropers = find_potential_impropers(bonds_with_names, donor_indices, metal_name)
    impropers = filter_aromatic_impropers(filename, potential_impropers, metal_name)
    return impropers


def find_dihedrals_and_values(bonds_with_names, metal_name, unique_ligands_pattern, starting_index, indecies, charge):
    print("Calculating improper dihedrals")

    impropers = find_impropers_connected_to_metal(bonds_with_names, metal_name)

    for unique_ligand in list(set(unique_ligands_pattern)):
        new_dihedrals = {}

        impropers_to_calculate = []

        same_ligands = np.where(unique_ligand==np.array(unique_ligands_pattern))[0]
        print(same_ligands)

        ligand_indecies = indecies[1:][same_ligands[0]]
        first_start_index = starting_index[1]

        for improper in impropers:
            if sum([1 for a in ligand_indecies if a in improper])==3:
                impropers_to_calculate.append(improper)
        print("impropers to calculate", impropers_to_calculate)

        # calculateing
        values_of_checked_impropers = []
        for improper in impropers_to_calculate:
            value = improper_value_calculation(improper, filename='bonded/site0_optimised_orca.xyz', charge=charge) #TODO we need to pass somethowe the name of optimised file

            values_of_checked_impropers.append(value)
            new_dihedrals[tuple(improper)] = tuple(value)
            new_dihedrals[tuple([improper[0],improper[2],improper[1],improper[3]])] = tuple(value)

        if len(same_ligands)>1:
            for same_ligand in same_ligands[1:]:
                print(same_ligand, starting_index[1:][same_ligand]-first_start_index)

                index_diffrence = starting_index[1:][same_ligand]-first_start_index
                proposed_torsions = np.array(impropers_to_calculate)+ np.array([index_diffrence, index_diffrence, index_diffrence, 0])
                #print(proposed_torsions)

                # check if found the same torsions:
                assert sum(1 for proposed_torsion in proposed_torsions if list(proposed_torsion) in impropers) == len(proposed_torsions)
                for idx_tor, proposed_torsion in enumerate(proposed_torsions):
                    new_dihedrals[tuple(proposed_torsion)] = tuple(values_of_checked_impropers[idx_tor])
                    new_dihedrals[tuple([proposed_torsion[0], proposed_torsion[2], proposed_torsion[1], proposed_torsion[3]])] = tuple(values_of_checked_impropers[idx_tor])

    return new_dihedrals

