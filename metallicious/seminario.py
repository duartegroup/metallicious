import os
import numpy as np

import MDAnalysis
from metallicious.mapping import map_two_structures
from metallicious.mod_seminario import modified_seminario_method
from metallicious.improper_torsion import find_impropers_and_values
from metallicious.log import logger

def check_if_orca_available():
    '''
    Checks if ORCA/autodE is available

    :return: (bool) True if ORCA/autodE is available
    '''
    try:
        import autode as ade
    except:
        raise ImportError("No autode found")

    method = ade.methods.ORCA()

    if method.is_available is False:
        raise NameError("For parametrization of templates, QM software ORCA is required")

def remove_non_metal_donor_bonds(bonds, metal_name, donors=['N', 'O', 'S']):
    '''
    Removes bonds which are not connected to donors (by default N, O, S). This is to remove bonds such as Metal-C or Metal-H

    :param bonds: (list) list of pairs of atoms which form bond (names)
    :param metal_name: (string) name of metal
    :param donors: (list) list of elements which metal can form a bond
    :return: (list) new list of pairs of atoms which form bond
    '''
    metal_name = metal_name.title()
    new_bonds = {}
    for bond in bonds:
        if metal_name.title() in bond[1]:
            if (bond[1][0] == metal_name and bond[1][1] in donors) or (
                    bond[1][1] == metal_name and bond[1][0] in donors):
                new_bonds[bond[0]] = bonds[bond]
            else:
                logger.info(
                    f"Seminario method detect bond {bond[1]}, but non-metal atom is not a donor (donors={donors:}), so removing it")

        else:
            new_bonds[bond[0]] = bonds[bond]
    return new_bonds


def symmetrize_bonds_angles(bonds, metal_name, filename_opt, indecies, unique_ligands_pattern):
    '''
    This procedure finds the bonds and angles which are equivalent in the graph representation

    :param bonds: (dict) dictionary including bonds parametrized by seminario
    :param metal_name: (str) metal name
    :param filename_opt: (str) path to optimised structure
    :param indecies: (list(list(int))) lists containing atom indices of different ligands
    :param unique_ligands_pattern: (list(int)) list showing which ligands are equivalent (e.g.,[0,0,1] -> first 2 are the same type of ligand). The bonded paramters of the same ligands are symmetrised
    :return: (dict) new dictionary including bond parameters
    '''

    site = MDAnalysis.Universe(filename_opt)
    ligand_indecies = indecies[1:]

    metal_index = 0
    for unique_ligand in list(set(unique_ligands_pattern)):
        unique_ligand_indecies = [a for a, b in enumerate(unique_ligands_pattern) if b == unique_ligand]

        mappings = []
        for unique_ligand_idx in unique_ligand_indecies:
            mapping, _ = map_two_structures(0, site.atoms[[0] + ligand_indecies[unique_ligand_idx]],
                                            site.atoms[[0] + ligand_indecies[unique_ligand_indecies[0]]], metal_name)
            mappings.append(mapping)

        bonds_reference = {}
        for bond in bonds:
            if len([1 for atom in bond if atom in [0] + ligand_indecies[unique_ligand_indecies[0]]]) == len(bond):
                bonds_reference[bond] = bonds[bond]

        for bond in bonds_reference:
            length = 0
            k = 0
            for mapping in mappings:
                mapping_idx = tuple([mapping[atom] for atom in bond])
                if mapping_idx in bonds:
                    mapped_bond = bonds[mapping_idx]
                elif mapping_idx[::-1] in bonds:
                    mapped_bond = bonds[mapping_idx[::-1]]

                length += mapped_bond[0]
                k += mapped_bond[1]

            length /= len(mappings)
            k /= len(mappings)

            for mapping in mappings:
                mapping_idx = tuple([mapping[atom] for atom in bond])

                if mapping_idx in bonds:
                    bonds[mapping_idx] = (np.round(length, 3), np.round(k, 3))
                elif mapping_idx[::-1] in bonds:
                    bonds[mapping_idx[::-1]] = (np.round(length, 3), np.round(k, 3))

    return bonds


def bond_remove_invalid_and_symmetrize(bonds_with_names, metal_name, filename_opt, indecies,
                                       unique_ligands_pattern, donors=["N", "O", "S"]):
    '''
    It starts procedures of removing not valid bonds (metal-non donor bonds, i.g., to prevent Metal-hydrogen).

    It normalizes values of the bonds which are the same

    :param bonds_with_names: (dict) dictionary including bonds
    :param metal_name: (str) name of the metal
    :param filename_opt: (str) filename of optimised structure
    :param starting_index: (list(int)) indices of first atoms of diffrent the residues
    :param indecies: (list(list)) list of indecies of diffrent residues
    :param unique_ligands_pattern: (list(int)) list showing which ligands are equivalent (e.g.,[0,0,1] -> first 2 are the same type of ligand). The bonded paramters of the same ligands are symmetrised
    :param donors: (list(str)) names of the donor atoms, with which the metal can form bond
    :return: (dict) dictionary which includes bonds
    '''

    bonds = remove_non_metal_donor_bonds(bonds_with_names, metal_name, donors=donors)
    bonds = symmetrize_bonds_angles(bonds, metal_name, filename_opt, indecies, unique_ligands_pattern)
    return bonds


def remove_non_metal_donor_angles(angles, metal_name, donors=['N', 'S', 'O']):
    '''
    It removes angle parameters which does not include bond metal-donor (for example to remove metal-hydrogen,
    or metal-carbon involving angles)

    :param angles: (dict) list of angle parameters
    :param metal_name: (str) name of the metal
    :param donors: (list(str)) list of names of donor atoms
    :return: (dict) new list of angle parameters with invalid angles removed
    '''
    metal_name = metal_name.title()
    new_angles = {}
    for angle in angles:
        if metal_name in angle[1]:
            metal_on_side_condition = (
                    (angle[1][0] == metal_name or angle[1][2] == metal_name) and (angle[1][1] in donors))
            metal_in_middle_condition = (
                    angle[1][1] == metal_name and (angle[1][0] in donors) and (angle[1][2] in donors))

            if metal_on_side_condition or metal_in_middle_condition:
                new_angles[angle[0]] = angles[angle]
            else:
                logger.info(
                    f"Not valid angle: Seminario method detected angle, but the atoms do not belong to donor list ({donors:}); angle:{angle:}")
        else:
            new_angles[angle[0]] = angles[angle]
    return new_angles


def angle_remove_invalid_and_symmetrize(angles_with_names, metal_name, filename_opt, indecies,
                                        unique_ligands_pattern,
                                        donors=["N", "O", "S"]):
    '''
    Removes angles from angles_with_names which include invalid metal-ligand bonded parameters, that is the one which
    are NOT formed with metal-donors atoms. In particular this is useful to remove weak bonds between metal and carbon,
    and metal and hydrogen.

    Then it normalizes the angles which are the same in different ligands.


    :param angles_with_names: (dict) dictionary of the angles
    :param metal_name: (str) name of the metal
    :param filename_opt: (str) filename of the optimized structure
    :param starting_index:
    :param indecies:
    :param unique_ligands_pattern: (list(int)) list showing which ligands are equivalent (e.g.,[0,0,1] -> first 2 are the same type of ligand). The bonded paramters of the same ligands are symmetrised
    :param donors:
    :return: (dict) new dictionary of the angles
    '''
    angles = remove_non_metal_donor_angles(angles_with_names, metal_name, donors=donors)
    angles = symmetrize_bonds_angles(angles, metal_name, filename_opt, indecies, unique_ligands_pattern)
    return angles


def extend_angle_to_dihedral(angle, bonds):
    '''
    Create a possible dihedral by extending angle by single bond

    :param angle: (list) indecies of atoms in angle
    :param bonds: (list) list of bonds
    :return: (list) dihedral
    '''
    dihedrals = []
    last_atom = angle[-1]
    temp = list(set(np.concatenate([[bond[0], bond[1]] for bond in bonds if last_atom in bond])))
    temp.remove(last_atom)
    for a in temp:
        if a not in angle:
            dihedrals.append(angle + [a])

    last_atom = angle[0]
    temp = list(set(np.concatenate([[bond[0], bond[1]] for bond in bonds if last_atom in bond])))
    temp.remove(last_atom)
    for a in temp:
        if a not in angle:
            dihedrals.append([a] + angle)
    return dihedrals


def generate_all_dihedrals(angles, bonds, metal_index=0):
    '''
    Create all possible dihedral by extending all angles by single bond

    :param angles: (list) list of angles
    :param bonds: (list) list of bonds
    :return: (list(list)) list of dihedrals
    '''
    dihedrals = []
    for angle_indexes in angles:
        if metal_index in angle_indexes:
            extended_dihedrals = extend_angle_to_dihedral(list(angle_indexes), bonds)
            dihedrals += extended_dihedrals

    # remove double counted dihedrals:
    new_dihedrals = []
    for a, dihedral in enumerate(dihedrals):
        if dihedral not in new_dihedrals and dihedral[::-1] not in new_dihedrals:
            new_dihedrals.append(dihedral)
    return new_dihedrals


def create_dummy_dihedrals(angles, bonds, metal_index=0):
    '''
    Generates dictionary of dihedrals which includes dihedrals generated by extending angles by single bond. The force
    constant and value of these are dihedrals are set to 0.

    :param angles: (list) list of angles
    :param bonds: (list) list of bonds
    :param metal_index: (int) metal index
    :return: (list(list)) list of dihedrals

    :return:
    '''
    dihedrals_indexes = generate_all_dihedrals(angles, bonds, metal_index)
    dihedrals = {}

    for dihedral in dihedrals_indexes:
        dihedrals[tuple(dihedral)] = (0, 0)
    return dihedrals


def create_pair_exclusions(dihedrals, angles):
    """
    Finds atoms interacting by 1-4 interactions, which will be put to pair exclusions
    :param dihedrals: (list(int,int,int,int)) list of indices quartets describing which atoms form dihedral
    :param angles: (list(int,int,int)) list of indices triads describing angles
    :return: list(int,int) list of indices pairs
    """
    pairs = []
    for dihedral in dihedrals:
        not_2_3_body = True

        for angle in angles:
            if dihedral[0] in angle and dihedral[-1] in angle:
                not_2_3_body = False

        if not_2_3_body:
            pairs.append((dihedral[0], dihedral[-1]))

    return pairs


def generate_angles_from_bonds(bond_list):
    """
    Finds angles by combining two bonds connecting the same atom
    :param bond_list: (list(int,int)) list of indices pairs describing bonds
    :return: (list(int,int,int)) list of indices triads describing angles
    """

    angle_list = []
    for bond in bond_list:
        angle_list += extend_angle_to_dihedral(list(bond), bond_list)
    return angle_list


def frequencies(filename, charge=0, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], mult=1):
    '''
    Runs frequency calculations using autodE

    :param filename: (str) name of the coordination file
    :param charge: (int) charge for QM calculations
    :param keywords: (list(str)) keywords for QM calculations
    :param mult: (int) multiplicity
    :return: (str, array(float)),list(str), array(array(float)), list((int,int)), list((int,int))
    name of optimised file, array containing position, list of names, hessian matrix, bond list (pairs of atoms), and
    angle list (trides of atoms)
    '''

    import autode as ade
    method = ade.methods.ORCA()

    if method.is_available is False:
        raise NameError("For parametrization of templates, QM software ORCA is required")

    site = ade.Molecule(filename, charge=charge, mult=mult)
    site.optimise(method=method, keywords=keywords)

    if site.imaginary_frequencies is not None:
        if len(site.imaginary_frequencies) > 0:
            # Sometimes it does not converge, then we use tighter criterium
            if 'TIGHTOPT' in [keyword.upper() for keyword in keywords] or 'OPT' in [keyword.upper() for keyword in
                                                                                    keywords]:
                new_keywords = []
                for keyword in keywords:
                    if keyword.upper() == 'TIGHTOPT' or keyword.upper() == 'OPT':
                        new_keywords.append("tightOPT")
                    else:
                        new_keywords.append(keyword)

                rerun_hessian = "\n%geom\nCalc_Hess true\nRecalc_Hess 10\nend"
                site.optimise(method=method, keywords=new_keywords + [rerun_hessian])
        else:
            raise RuntimeError("Optimization fails")

    names = [atom.atomic_symbol for atom in site.atoms]
    opt_filename = f"{site.name:s}_optimised.xyz"
    site.print_xyz_file(filename=opt_filename)

    # we use .to(ade.units.J_per_ang_sq) twice due to changes from autode 1.3.2 to 1.4.0
    # this is for compability with autode 1.3.2
    site.hessian.to(ade.units.J_per_ang_sq)
    # and here it is for compability with autod 1.4.0
    hessian = np.array(site.hessian.to(ade.units.J_per_ang_sq))
    hessian *= 6.022e23 / 4.184e3  # we convert tto kcal/mol/angs^2

    bond_list = list(site.graph.edges)
    angle_list = generate_angles_from_bonds(bond_list)

    return opt_filename, np.array(site.coordinates), names, hessian, bond_list, angle_list


def simple_seminario(filename, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], charge=0, mult=1,
                     vibrational_scaling=None):
    '''
    Starts frequency calculations and the uses Mod. Seminario method to get bonds and angles

    :param filename: (str) input structure file
    :param keywords: (list(str)) keywords for QM calculations
    :param charge: (int) charge for QM calculations
    :param mult: (int) multiplicity for charge calculations
    :param vibrational_scaling: (float) vibrational scaling, see  https://cccbdb.nist.gov/vsfx.asp
    :return:  (dict(), list(), str) two dictionaries, bonds and angles, with indices of atoms and their names as key,
    and bond lenght and force constant as values; name of the optimised file
    '''

    here = os.getcwd()
    os.system("mkdir bonded")
    os.chdir("bonded")

    opt_filename, coords, atom_names, hessian, bond_list, angle_list = frequencies("../" + filename, charge,
                                                                                   keywords=keywords, mult=mult)

    # Scalling factor taken from https://cccbdb.nist.gov/vsfx.asp and assumed that for basis set of double zeta and higher the scalling does not change
    if vibrational_scaling is not None:
        vibrational_scaling = vibrational_scaling
    elif 'PBE0' in keywords:
        vibrational_scaling = 0.96
    elif 'wB97X-D3' in keywords:
        vibrational_scaling = 0.955
    else:
        vibrational_scaling = 1

    bonds_with_names, angles_with_names = modified_seminario_method(hessian, coords, atom_names, bond_list, angle_list,
                                                                    vibrational_scaling=vibrational_scaling)

    os.chdir(here)

    return bonds_with_names, angles_with_names, f"bonded/{opt_filename:s}"


def remove_atoms_from_bonded(bonds, angles, impropers, atoms_to_remove):
    '''
    Removes bonds, angles, and impropers which include specified atoms (atoms_to_remove). This procedure is used to
    remove hydrogen atoms needed for the QM calculations, but removed to create template.

    :param bonds: (dict) dictionary of the bonds
    :param angles: (dict) dictionary of the angles
    :param impropers: (dict) dictionary of the impropers
    :param atoms_to_remove: (list(int)) list of indices of atoms
    :return: (dict, dict, dict) new dictionaries of bonds, angles and impropers without parameters for specified atoms
    '''
    new_params = []
    for bonded in [bonds, angles, impropers]:
        new_bonded = {}
        for param in bonded:
            remove = False
            for atom_id in param:
                if atom_id in atoms_to_remove:
                    remove = True

            if not remove:
                new_bonded[param] = bonded[param]

        new_params.append(new_bonded)

    new_bonds, new_angles, new_impropers = new_params
    return new_bonds, new_angles, new_impropers


def remove_indices_of_removed_atoms(bonds, angles, impropers, atoms_to_remove):
    '''
    Renumbers the atom indices in bonded parameters that thy are in a sequence without gaps

    :param bonds: (dict) dictionary of the bonds
    :param angles: (dict) dictionary of the angles
    :param impropers: (dict) dictionary of the impropers
    :param atoms_to_remove: (list(int)) list of indices of atoms
    :return: (dict, dict, dict) new dictionaries of bonds, angles and impropers with reordered
    '''
    atoms_to_remove = sorted(atoms_to_remove)
    new_params = []
    for bonded in [bonds, angles, impropers]:
        new_bonded = {}
        for param in bonded:
            indices = tuple([idx - sum(np.array(atoms_to_remove) < idx) for idx in param])
            new_bonded[indices] = bonded[param]
        new_params.append(new_bonded)
    new_bonds, new_angles, new_impropers = new_params

    return new_bonds, new_angles, new_impropers


def single_seminario(filename, metal_charge, metal_name, starting_index, indecies, unique_ligands_pattern,
                     keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], mult=1, improper_metal=False,
                     donors=['N', 'S', 'O'], atoms_to_remove=None, vibrational_scaling=None):
    """
    Runs seminario method for the structure inputed in filename

    :param filename: (str) input structure for which bonds and angles will be found
    :param metal_charge: (int) charge for QM calculations
    :param metal_name: (str) name of the metal
    :param starting_index: list(int) indices where the new ligands start
    :param indecies: (list(list(int))) lists containing atom indices of different ligands
    :param unique_ligands_pattern: (list(int)) list showing which ligands are equivalent (e.g.,[0,0,1] -> first 2 are the same type of ligand). The bonded paramters of the same ligands are symmetrised
    :param keywords: (list(str)) keywords for level of theory
    :param mult: (int) multiplitiy of the metal
    :param improper_metal: (bool) if True, parametrize the improper dihedral involving metal
    :param donors: (list(str)) atoms with which metal can interact
    :param atoms_to_remove: (list) list of hydrogen atoms to remove in order to create the template from saturated template
    :param vibrational_scaling: (float) vibrational scaling for the QM frequencies
    :return:
    """

    # bonds_with_names, angles_with_names, dummy_dihedrals, filename_opt = simple_seminario(filename, keywords=keywords, charge=metal_charge, mult=mult, vibrational_scaling=vibrational_scaling)
    bonds_with_names, angles_with_names, filename_opt = simple_seminario(filename, keywords=keywords,
                                                                         charge=metal_charge,
                                                                         mult=mult,
                                                                         vibrational_scaling=vibrational_scaling)

    bonds = bond_remove_invalid_and_symmetrize(bonds_with_names, metal_name, filename_opt, indecies,
                                               unique_ligands_pattern,
                                               donors=donors)

    angles = angle_remove_invalid_and_symmetrize(angles_with_names, metal_name, filename_opt, indecies,
                                                 unique_ligands_pattern,
                                                 donors=donors)

    if improper_metal:
        impropers = find_impropers_and_values(bonds, metal_name, unique_ligands_pattern, starting_index, indecies,
                                              charge=metal_charge, mult=mult, filename=filename_opt)
    else:
        impropers = {}

    if len(atoms_to_remove) is not None:
        bonds, angles, impropers = remove_atoms_from_bonded(bonds, angles, impropers, atoms_to_remove)
        bonds, angles, impropers = remove_indices_of_removed_atoms(bonds, angles, impropers, atoms_to_remove)

    dummy_dihedrals = create_dummy_dihedrals(angles, bonds, metal_index=0)

    dihedrals = {}  # dihedrals are not implemented

    pairs = create_pair_exclusions(dummy_dihedrals, angles)

    return bonds, angles, dihedrals, impropers, pairs, filename_opt
