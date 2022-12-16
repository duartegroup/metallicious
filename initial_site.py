from antechamber_interface import antechamber
import parmed as pmd


try:
    from cgbind2pmd.data import name2mass
    from cgbind2pmd.data import name_to_atomic_number
    from cgbind2pmd.data import vdw_data
except:
    from data import name2mass
    from data import name_to_atomic_number
    from data import vdw_data

def create_metal_topol(metal_name, metal_charge, vdw_data_name): #TODO

    print(vdw_data)
    data = vdw_data[vdw_data_name]

    if metal_name.title() in data:
        eps = data[metal_name][0] # TODO check units
        radius = data[metal_name][1]

        mass = name2mass[metal_name] # TODO refactor names
        atomic_number = name_to_atomic_number[metal_name]

        new_atom = pmd.topologyobjects.Atom(atomic_number=atomic_number, type=metal_name, name=metal_name,
                                            rmin=radius, epsilon=eps, mass=mass, charge=metal_charge)
        struct = pmd.structure.Structure()
        struct.add_atom(new_atom, metal_name, metal_name)

        # I don't know how to convert struct to GromacsTopol, easiest is to save and load it ...
        struct.save(f"{metal_name:s}.top", overwrite=True)
        metal_topol = pmd.load_file(f"{metal_name:s}.top")
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
                antechamber(f"ligand{n_site:d}_{idx:}.pdb", 0, f"ligand{n_site:d}_{idx:}.itp") #TODO n_site sometimes called idx_site
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



