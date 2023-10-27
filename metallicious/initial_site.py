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
        raise Exception(f"Metal {metal_name.title():}{metal_charge:} not in the library")

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
        metal_topol.defaults.fudgeLJ = 0.5
        metal_topol.defaults.fudgeQQ = 0.833333

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


