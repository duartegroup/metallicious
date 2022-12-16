import argparse

from extract_metal_site import extract_metal_structure
from seminario import multi_seminario
from charges import calculate_charges
from initial_site import create_initial_topol

import parmed as pmd
import re
import numpy as np

import MDAnalysis

def strip_numbers_from_atom_name(atom_name):
    return re.match("([a-zA-Z]+)", atom_name).group(0)

def copy_bonds(topol_new, bonds, metal_name):
    '''
    Copies bonds into topology

    :param topol_new:
    :param bonds:
    :return:
    '''
    orginal_bonds = []

    for bond in topol_new.bonds:
        orginal_bonds.append((bond.atom1.idx, bond.atom2.idx))


    for bond in bonds:
        find = [a for a, bond_o in enumerate(orginal_bonds) if (bond_o == bond or bond_o == bond[::-1])]
        print(bond, find)
        if len(find) == 1:
            bond_new = topol_new.bonds[find[0]]

            type_to_assign = pmd.topologyobjects.BondType(bonds[bond][1], bonds[bond][0], list=topol_new.bond_types)
            # if type_to_assign not in self.topol_new.bond_types:
            topol_new.bond_types.append(type_to_assign)

            # bond_new.type = deepcopy(bond_fp.type)
            bond_new.type = type_to_assign  # bond_fp.type
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[bond[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                    topol_new.atoms[bond[1]].name).title() == metal_name.title():
                type_to_assign = pmd.topologyobjects.BondType(bonds[bond][1], bonds[bond][0], list=topol_new.bond_types)
                topol_new.bond_types.append(type_to_assign)
                atom1 = topol_new.atoms[bond[0]]
                atom2 = topol_new.atoms[bond[1]]
                topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))
            else:
                print("this is wierd bond", topol_new.atoms[bond[0]].name, topol_new.atoms[bond[1]].name, bond,
                      bonds[bond])
    return topol_new

def copy_angles(topol_new, angles, metal_name):
    '''
    Copy angles into topology

    :param topol_new:
    :param angles:
    :return:
    '''
    orginal_angles = []

    for angle in topol_new.angles:
        orginal_angles.append((angle.atom1.idx, angle.atom2.idx, angle.atom3.idx))

    for angle in angles:
        find = [a for a, angle_o in enumerate(orginal_angles) if (angle_o == angle or angle_o == angle[::-1])]
        if len(find) == 1:
            angle_new = topol_new.angles[find[0]]
            type_to_assign = pmd.topologyobjects.AngleType(angles[angle][1], angles[angle][0],
                                                           list=topol_new.angle_types)

            topol_new.angle_types.append(type_to_assign)
            angle_new.type = type_to_assign
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[angle[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                    topol_new.atoms[angle[1]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                    topol_new.atoms[angle[2]].name).title() == metal_name.title():

                type_to_assign = pmd.topologyobjects.AngleType(angles[angle][1], angles[angle][0],
                                                               list=topol_new.angle_types)

                topol_new.angle_types.append(type_to_assign)

                atom1 = topol_new.atoms[angle[0]]
                atom2 = topol_new.atoms[angle[1]]
                atom3 = topol_new.atoms[angle[2]]

                topol_new.angles.append(pmd.topologyobjects.Angle(atom1, atom2, atom3, type=type_to_assign))
            else:
                print('nope', angle)

    return topol_new


def copy_dihedrals(topol_new, dihedrals, metal_name):
    '''
    Copy angles into topology

    :param topol_new:
    :param angles:
    :return:
    '''
    orginal_dihedrals = []

    for dihedral in topol_new.dihedrals:
        orginal_dihedrals.append((dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx))

    for dihedral in dihedrals:
        find = [a for a, dihedral_o in enumerate(orginal_dihedrals) if
                (dihedral_o == dihedral or dihedral_o == dihedral[::-1])]
        if len(find) == 1:
            print("It found dihedral in the topology, but it should not, as they should be not present, doing noting")
            print(dihedral)  # TODO
            raise Error
            # This procedures are not done in reality
            # dihedral_new = topol_new.dihedrals[find[0]]
            # type_to_assign = pmd.topologyobjects.DihedralType(dihedrals[dihedral][1], dihedrals[dihedral][0],list=topol_new.angle_types)
            # topol_new.angle_types.append(type_to_assign)
            # dihedral_new.type = type_to_assign
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[dihedral[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[1]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[2]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[3]].name).title() == metal_name.title():

                # we do use periodicy=2 and use defualt values for scee=1.2 and scnb=2 (screening of electrostatics and noboned)
                type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedrals[dihedral][1], per=2,
                                                                  phase=dihedrals[dihedral][0],
                                                                  list=topol_new.dihedral_types)

                topol_new.dihedral_types.append(type_to_assign)

                atom1 = topol_new.atoms[dihedral[0]]
                atom2 = topol_new.atoms[dihedral[1]]
                atom3 = topol_new.atoms[dihedral[2]]
                atom4 = topol_new.atoms[dihedral[3]]
                topol_new.dihedrals.append(
                    pmd.topologyobjects.Dihedral(atom1, atom2, atom3, atom4, type=type_to_assign))
            else:
                print('nope', angle)

    return topol_new


def parametrize(filename, metal_name, metal_charge, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], vdw_type="uff", mult=1):
    '''
    Combines all other scripts to extract, and paramterize metal site(s)
    :param filename:
    :param metal_name:
    :param metal_charge:
    :return:
    '''
    print("[ ] Extracting the structure")
    extract_metal_structure(filename, metal_name, output="site")
    print("[ ] Create initial topology")
    topols = create_initial_topol(metal_name, metal_charge, vdw_data_name=vdw_type)
    n_sites = len(topols)

    topols = []
    for n_site in range(n_sites):
        topols.append(pmd.load_file(f"topol{n_site:d}.top"))
    print(topols)
    print(topols[0])

    print("[ ] Peforming seminario calculations")
    bonds, angles, dummy_dihedrals = multi_seminario(metal_charge, metal_name, keywords=keywords, mult=mult)
    print("[ ] Seminario finished, now coping")
    print(bonds)
    print("------------------")
    for a, topol in enumerate(topols):
        topols[a] = copy_bonds(topols[a], bonds[a], metal_name)
        topols[a] = copy_angles(topols[a], angles[a], metal_name)
        topols[a] = copy_dihedrals(topols[a], dummy_dihedrals[a], metal_name)

    print("Calcualting charges")
    partial_charges = calculate_charges(metal_charge, metal_name,vdw_data_name=vdw_type, mult=mult)

    print("Copying charges")
    for n_site, topol in enumerate(topols):
        for idx, atom in enumerate(topol.atoms):
            atom.charge = partial_charges[n_site][idx]
    print("Charges calculated !")

    File = open("INFO.dat")
    text = File.read()
    File.close()
    for line in text.splitlines():
        if "extra_atoms:" in line:
            extra_atoms = list(map(int, line[12:].split(',')))

    for a, topol in enumerate(topols):
        topol.strip(f"@{','.join(list(map(str, np.array(extra_atoms) + 2))):s}")
        topol.write(f"new_topol{a:d}.top")

        new_cage = MDAnalysis.Universe(f"site{a:d}.pdb")
        new_cage.select_atoms(f"not index {' '.join(list(map(str, np.array(extra_atoms) + 1))):s}").write(f"template{a:d}.pdb")


    print("Cut the hydrogens")


    print("Program finishess successfully!")

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-o", help="Clean structure")
    parser.add_argument("-metal_name", help="Name of the metal")
    parser.add_argument("-metal_charge", help="Charge of the metal")
    parser.add_argument("-keywords", help="keywords for QM", nargs='+')
    parser.add_argument("-vdw_type", default='uff', help="Type of parameters for VdW (available: uff, merz-tip3p, merz-opc3, merz-spc/e, merz-tip3p-fb, merz-opc, merz-tip4p-fb, merz-tip4-ew, zhang-tip3p, zhang-opc3, zhang-spc/e, zhang-spc/eb, zhang-tip3p-fb, zhang-opc, zhang-tip4p/2005, zhang-tip4p-d, zhang-tip4p-fb, zhang-tip4p-ew")
    parser.add_argument("-mult", default=1, help="multiplicity")
    return parser.parse_args()




if __name__ == '__main__':
    args = get_args()
    #print("AAAA", args.keywords)
    parametrize(args.f, args.metal_name, int(args.metal_charge), args.keywords, vdw_type=args.vdw_type, mult=int(args.mult))


