import parmed as pmd
import os

# try:
from metallicious.data import name2mass
from metallicious.data import name_to_atomic_number
from metallicious.data import vdw_data
from metallicious.antechamber_interface import antechamber
# except:
#     from data import name2mass
#     from data import name_to_atomic_number
#     from data import vdw_data
#     from antechamber_interface import antechamber


def create_metal_topol(metal_name, metal_charge, vdw_data_name):
    data = vdw_data[vdw_data_name]

    if f"{metal_name.title():}{metal_charge:}" in data:
        name = f"{metal_name.title():}{metal_charge:}"
    elif metal_name in data:
        name = metal_name
    else:
        print(f"Metal {metal_name.title():}{metal_charge:} not in the library")
        raise

    if name in data:
        eps = data[name][0]
        radius = data[name][1]

        mass = name2mass[metal_name]
        atomic_number = name_to_atomic_number[metal_name]

        new_atom = pmd.topologyobjects.Atom(atomic_number=atomic_number, type=metal_name, name=metal_name,
                                            rmin=radius, epsilon=eps, mass=mass, charge=metal_charge, number=0)

        struct = pmd.structure.Structure()
        struct.add_atom(new_atom, metal_name, metal_name)

        # I don't know how to convert struct to GromacsTopol, easiest is to save and load it ...
        struct.save(f"{name:s}.top", overwrite=True)
        metal_topol = pmd.load_file(f"{name:s}.top")
        os.remove(f"{name:s}.top")

        ''' I have several attemps how to make it without file but did not manage to figure out:
        metal_topol = pmd.gromacs.GromacsTopologyFile()
        metal_topol.add_atom(new_atom, metal_name, metal_name)
        #atomtype = pmd.topologyobjects.AtomType("Pd", None, 106)

        residue = pmd.topologyobjects.Residue(metal_name, number=0)
        residue.add_atom(new_atom)
        residue_list = pmd.topologyobjects.ResidueList()
        residue_list.append(residue)
        metal_topol.residues = residue_list
        '''

        return metal_topol


''' TODO me trying to find out how to do this thing
prop = [p for p in dir(pmd.topologyobjects.Atom) if isinstance(getattr(pmd.topologyobjects.Atom,p),property)]
prop = [p for p in dir(pmd.gromacs.gromacstop.GromacsTopologyFile) if isinstance(getattr(pmd.gromacs.gromacstop.GromacsTopologyFile,p),property)]

atom1=metal
atom2=metal2

for pro in prop:
    if hasattr(atom1, pro) and hasattr(atom2, pro):
        print("SAME:", getattr(atom1,pro), getattr(atom2,pro))
    elif hasattr(atom1, pro):
        print("atom1 only:", getattr(atom1, pro))
    elif hasattr(atom2, pro):
        print("atom2 only:", getattr(atom2, pro))
    else:
        print("nope", pro)

'''
def create_initial_topol(metal_name, metal_charge, vdw_data_name):
    '''
    This creates initial (covalent) topology consistent of the metal ions and residues

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
                antechamber(f"ligand_{idx:}.pdb", f"ligand_{idx:}.itp")
                unique_ligands_topologies.append(pmd.load_file(f"ligand{n_site:d}_{idx:}.itp"))

            metal_topol = create_metal_topol(metal_name, metal_charge, vdw_data_name)

            n_metals = 1

            topol = metal_topol


            for pattern in unique_ligands_pattern:
                print(pattern)
                topol += unique_ligands_topologies[pattern]

            topol.write(f"topol.top", [list(range(n_metals + len(unique_ligands_pattern)))])
            topols.append(topol)

            n_site+=1
    return topols

def create_initial_topol2(metal_name, metal_charge, unique_ligands_pattern, vdw_data_name):
    unique_ligands_topologies = []
    unique_ligands = list(set(unique_ligands_pattern))

    for idx, unique_ligand in enumerate(unique_ligands):
        antechamber(f"ligand_{idx:}.pdb", f"ligand_{idx:}.itp")
        unique_ligands_topologies.append(pmd.load_file(f"ligand_{idx:}.itp"))

    metal_topol = create_metal_topol(metal_name, metal_charge, vdw_data_name)

    n_metals = 1
    topol = metal_topol

    for pattern in unique_ligands_pattern:
        topol += unique_ligands_topologies[pattern]

    topol.write(f"topol.top", [list(range(n_metals + len(unique_ligands_pattern)))])
    return topol

def create_initial_topol_from_larger_topol(metal_name, metal_charge, unique_ligands_pattern, vdw_data_name):
    unique_ligands_topologies = []
    unique_ligands = list(set(unique_ligands_pattern))

    for idx, unique_ligand in enumerate(unique_ligands):
        antechamber(f"ligand_{idx:}.pdb", f"ligand_{idx:}.itp")
        unique_ligands_topologies.append(pmd.load_file(f"ligand_{idx:}.itp"))

    metal_topol = create_metal_topol(metal_name, metal_charge, vdw_data_name)

    n_metals = 1
    topol = metal_topol

    for pattern in unique_ligands_pattern:
        topol += unique_ligands_topologies[pattern]

    topol.write(f"topol.top", [list(range(n_metals + len(unique_ligands_pattern)))])
    return topol


import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-metal_charge", help="Clean structure")
    parser.add_argument("-metal_name", help="Clean structure")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    create_initial_topol()



