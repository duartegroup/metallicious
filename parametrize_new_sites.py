import argparse
# TODO this should be class not INFO file
from extract_metal_site import extract_metal_structure
from seminario import multi_seminario
from charges import calculate_charges
from initial_site import create_initial_topol
from copy_topology_params import copy_bonds, copy_angles, copy_dihedrals

import parmed as pmd
import numpy as np
import MDAnalysis
import os


def parametrize(filename, metal_name, metal_charge, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], vdw_type="uff", mult=1, improper_metal=False, donors = ['N', 'S', 'O']):
    '''
    Combines all other scripts to extract, and paramterize metal site(s)
    :param filename:
    :param metal_name:
    :param metal_charge:
    :return:
    '''

    if os.path.isfile("INFO.dat"):
        os.remove('INFO.dat')


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
    bonds, angles, dummy_dihedrals = multi_seminario(metal_charge, metal_name, keywords=keywords, mult=mult, improper_metal=improper_metal, donors=donors)
    print("\t[ ] Seminario finished, now coping")
    for a, topol in enumerate(topols):
        topols[a] = copy_bonds(topols[a], bonds[a], metal_name)
        topols[a] = copy_angles(topols[a], angles[a], metal_name)
        topols[a] = copy_dihedrals(topols[a], dummy_dihedrals[a], metal_name)
    print("\t[+] Copied!")


    print("[ ] Calcualting charges")
    partial_charges = calculate_charges(metal_charge, metal_name,vdw_data_name=vdw_type, mult=mult)

    print("\t[ ] Copying charges")
    for n_site, topol in enumerate(topols):
        for idx, atom in enumerate(topol.atoms):
            atom.charge = partial_charges[n_site][idx]
    print("\t[+] Charges calculated !")


    print("[ ] Cuting extra atoms")
    File = open("INFO.dat") # TODO this here is wierd
    text = File.read()
    File.close()
    for line in text.splitlines():
        if "extra_atoms:" in line:
            extra_atoms = list(map(int, line[12:].split(',')))


    for a, topol in enumerate(topols):
        topol.write(f"old_new_topol{a:d}.top")

        topol.strip(f"@{','.join(list(map(str, np.array(extra_atoms) + 1))):s}")
        topol.write(f"new_topol{a:d}.top")
        new_cage = MDAnalysis.Universe(f"site{a:d}.pdb")
        new_cage.atoms.write(f"old_new_template{a:d}.pdb")

        print(f"We use this on site{a:d}.pdb:")
        print(f"not index {' '.join(list(map(str, np.array(extra_atoms)))):s}")
        new_cage.select_atoms(f"not index {' '.join(list(map(str, np.array(extra_atoms)))):s}").write(f"template{a:d}.pdb")

    print("Hydrogens removed!")


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
    parser.add_argument("-improper_metal", action='store_true', default=False, help="Calculate improper dihedral of the metal-aromatic (default:False)")
    parser.add_argument("-donors", nargs ='+', default=['N', 'S', 'O'],
                        help="Donors from the connected ligands, usually electronegative atom, such as N, S, O, but sometimes metal is connected to carbon", )
    return parser.parse_args()




if __name__ == '__main__':
    args = get_args()
    #print("AAAA", args.keywords)
    parametrize(args.f, args.metal_name, int(args.metal_charge), args.keywords, vdw_type=args.vdw_type, mult=int(args.mult), improper_metal=args.improper_metal, donors = args.donors)


