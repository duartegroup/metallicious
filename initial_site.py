import parmed as pmd


try:
    from cgbind2pmd.data import name2mass
    from cgbind2pmd.data import name_to_atomic_number
    from cgbind2pmd.data import vdw_data
    from cgbind2pmd.antechamber_interface import antechamber
except:
    from data import name2mass
    from data import name_to_atomic_number
    from data import vdw_data
    from antechamber_interface import antechamber


def create_metal_topol(metal_name, metal_charge, vdw_data_name):
    data = vdw_data[vdw_data_name]

    if metal_name.title() in data:
        eps = data[metal_name][0]
        radius = data[metal_name][1]

        mass = name2mass[metal_name]
        atomic_number = name_to_atomic_number[metal_name]

        new_atom = pmd.topologyobjects.Atom(atomic_number=atomic_number, type=metal_name, name=metal_name,
                                            rmin=radius, epsilon=eps, mass=mass, charge=metal_charge)
        # old -- to be removed when tested
        #struct = pmd.structure.Structure()
        #struct.add_atom(new_atom, metal_name, metal_name)

        # I don't know how to convert struct to GromacsTopol, easiest is to save and load it ...
        #struct.save(f"{metal_name:s}.top", overwrite=True)
        #metal_topol = pmd.load_file(f"{metal_name:s}.top")
        # end of old - above to be removed

        metal_topol = pmd.gromacs.GromacsTopologyFile()
        metal_topol.add_atom(new_atom, metal_name, metal_name)

        residue = pmd.topologyobjects.Residue(metal_name, number=0)
        residue.add_atom(new_atom)
        residue_list = pmd.topologyobjects.ResidueList()
        residue_list.append(residue)
        metal_topol.residues = residue_list

        return metal_topol
    else:
        print("Metal ", metal_name.title(), "not in the library")
        raise

def create_initial_topol(metal_name, metal_charge, vdw_data_name):
    '''
    This creates initial (noncovalentl) topology consistent of the metal ions and residues

    :param metal_name:
    :return:
    '''

    File = open("INFO.dat")
    text = File.read()
    File.close()

    unique_ligands_pattern = None

    all_sites_charges = []

    n_site = 0
    topols = []
    for line in text.splitlines():
        if "ligand_pattern:" in line:
            unique_ligands_pattern = list(map(int, line[15:].split(',')))
            unique_ligands = list(set(unique_ligands_pattern))

            unique_ligands_topologies = []

            for idx, unique_ligand in enumerate(unique_ligands):
                antechamber(f"ligand{n_site:d}_{idx:}.pdb", 0, f"ligand{n_site:d}_{idx:}.itp")
                unique_ligands_topologies.append(pmd.load_file(f"ligand{n_site:d}_{idx:}.itp"))

            metal_topol = create_metal_topol(metal_name, metal_charge, vdw_data_name)

            n_metals = 1

            topol = metal_topol


            for pattern in unique_ligands_pattern:
                print(pattern)
                topol += unique_ligands_topologies[pattern]

            topol.write(f"topol{n_site:d}.top", [list(range(n_metals + len(unique_ligands_pattern)))])
            topols.append(topol)

        n_site+=1
    return topols

import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-metal_charge", help="Clean structure")
    parser.add_argument("-metal_name", help="Clean structure")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    create_initial_topol()



