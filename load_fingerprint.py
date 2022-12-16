import os
import parmed as pmd
import MDAnalysis
import numpy as np


def load_fingerprint(name_of_binding_side, fingerprint_style):
    topol = pmd.load_file(f'{os.path.dirname(__file__):s}/library/{name_of_binding_side:s}.top')
    syst_fingerprint = MDAnalysis.Universe(f"{os.path.dirname(__file__):s}/library/{name_of_binding_side:s}.pdb")

    if fingerprint_style=='full':
        return topol, syst_fingerprint

    elif fingerprint_style=='trunc':
        metal_topology = topol.atoms[0]
        atoms_bound_to_metal_by_angle = list(set(
            np.concatenate([[angle.atom1.idx, angle.atom2.idx, angle.atom3.idx] for angle in metal_topology.angles])))
        atoms_bound_to_metal_by_bond = list(
            set(np.concatenate([[bond.atom1.idx, bond.atom2.idx] for bond in metal_topology.bonds])))

        atoms_to_strip = []
        for atom in topol.atoms:
            if atom.idx not in atoms_bound_to_metal_by_angle:
                atoms_to_strip.append(atom.idx)

        residual_charge = 0.0
        for atom in topol.atoms:
            if atom.idx not in atoms_bound_to_metal_by_bond:
                residual_charge += atom.charge
        residual_charge /= len(atoms_bound_to_metal_by_bond)

        for idx in atoms_bound_to_metal_by_angle:
            if idx in atoms_bound_to_metal_by_bond:
                topol.atoms[idx].charge += residual_charge
            else:
                topol.atoms[idx].charge = 0.0

        new_topol = topol.strip(f"@{','.join(list(map(str, np.array(atoms_to_strip) + 1))):s}")
        new_syst_fingerprint = syst_fingerprint.select_atoms(f"not index {' '.join(list(map(str, np.array(atoms_to_strip)))):s}")
        return new_topol, new_syst_fingerprint

    elif fingerprint_style=='min':
        metal_topology = topol.atoms[0]
        atoms_bound_to_metal_by_bond = list(
            set(np.concatenate([[bond.atom1.idx, bond.atom2.idx] for bond in metal_topology.bonds])))
        atoms_to_strip = []
        for atom in topol.atoms:
            if atom.idx not in atoms_bound_to_metal_by_bond:
                atoms_to_strip.append(atom.idx)

        for atom in topol.atoms:
            atom.charge = 0.0

        new_topol = topol.strip(f"@{','.join(list(map(str, np.array(atoms_to_strip) + 1))):s}")
        new_syst_fingerprint = syst_fingerprint.select_atoms(f"not index {' '.join(list(map(str, np.array(atoms_to_strip)))):s}")
        return new_topol, new_syst_fingerprint