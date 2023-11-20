import argparse
from subprocess import Popen, PIPE
import numpy as np
import os
import shutil
from metallicious.log import logger
from metallicious.utils import new_directory

import MDAnalysis
def connectivity_for_qcel(filename, metal_name):
    atoms = MDAnalysis.Universe(filename).atoms
    if metal_name is not None:
        types = atoms.types
        types[0] = metal_name
        atoms.types = types
    bonds = MDAnalysis.topology.guessers.guess_bonds(atoms, atoms.positions, vdwradii={metal_name: 0.1})
    new_bonds = [(bond[0], bond[1], 1) for bond in bonds]
    return new_bonds


def resp_orca(filename, charge=0, opt=True, metal_name=None, metal_radius=None, n_reorientations=1, mult=1,
              extra_atoms=None, keywords = ['D3BJ', 'PBE0', 'def2-SVP', 'keepdens']):
    '''
    Single point calculations of ESP and RESP analysis using autode/ORCA and psiresp
    This part of the script overwrites parts of the functionality of autode

    :param filename: (string) name of *xyz file
    :param charge: (int) charge of the whole structure
    :param opt: (boolean) optimise structure before single point
    :param metal_name: (string) name of the matal ion
    :param metal_radius: (float) VdW radius of the metal
    :param n_reorientations: (int) how many diffrent ortientation of ESP calculate
    :param mult: (int) multiplicity of the system
    :param extra_atoms: (list(int)) group of atoms which will be constrained to charge zero
    :param keywords: (list(str)) keywords for the autode
    :return:
    '''


    import autode as ade
    from autode.wrappers.ORCA import work_in_tmp_dir, logger
    import psiresp
    import qcelemental as qcel

    # Below functions are to overwrite some of the autode functionality
    # in particular we want to save the ESP files, densities etc.
    def old_execute(self, calc):
        @work_in_tmp_dir(filenames_to_copy=calc.input.filenames,
                         kept_file_exts=('.out', '.hess', '.xyz', '.inp', '.pc'))
        def execute_orca():
            # Run the calculations
            run_external(params=[calc.method.path, calc.input.filename],
                         output_filename=calc.output.filename)

        execute_orca()
        return None

    # @check_sufficient_memory
    def run_external(params, output_filename):
        """
        Standard method to run a EST calculation with subprocess writing the
        output to the calculation output filename

        ---------------------------------------------------------------------------
        Arguments:
            output_filename (str):

            params (list(str)): e.g. [/path/to/method, input-filename]
        """
        print("running new execute", params)
        with open(output_filename, 'w') as output_file:
            # /path/to/method input_filename > output_filename
            process = Popen(params, stdout=output_file, stderr=PIPE, close_fds=True)

            with process.stderr:
                for line in iter(process.stderr.readline, b''):
                    logger.warning('STDERR: %r', line.decode())
                    print(line.decode())

            poll = process.poll()
            if poll is None:
                print("Process is running 1")
            process.wait()

            poll = process.poll()
            if poll is None:
                print("Process is running 2")
        return None

    def new_execute(self, calc):
        @work_in_tmp_dir(filenames_to_copy=calc.input.filenames + ["grid.vpot.xyz"],
                         kept_file_exts=('.out', '.hess', '.xyz', '.inp', '.pc', '.input', '.densities', '.gbw'))
        def execute_orca():
            # Run the calculations
            run_external(params=[calc.method.path, calc.input.filename],
                         output_filename=calc.output.filename)
            # orca_vpot filename.gbw filename.scfp filename.vpot.xyz filename.vpot.out
            filename_core = calc.input.filename.replace(".inp", "")
            run_external(
                params=[calc.method.path + "_vpot", filename_core + ".gbw", filename_core + ".scfp", "grid.vpot.xyz",
                        filename_core + ".vpot.out"],
                output_filename=filename_core + "_vpot.out")
        execute_orca()
        return None

    method = ade.methods.ORCA()
    if method.is_available is False:
        raise NameError("For parametrization of templates, QM software ORCA is required")

    site = ade.Molecule(filename, charge=charge, mult=mult)
    molecule_qcel = qcel.models.Molecule.from_file(filename, molecular_charge=charge, connectivity=connectivity_for_qcel(filename, metal_name))
    molecule_psiresp = psiresp.Molecule(qcmol=molecule_qcel, charge=charge)
    logger.info(f"Formal charge of the molecule for template calculations: {charge:}")

    here = os.getcwd()

    # Becasue how autode works, it will try to load the files, but becasue we modified how autode is executed here
    # t might load incorrect grid
    # we need to remove folder with previous calculations if they exist:

    if os.path.isdir("resp"):
        shutil.rmtree("resp")
    new_directory("resp")
    os.chdir("resp")


    if opt:
        ade.wrappers.ORCA.ORCA.execute = old_execute
        site.optimise(method=method)

    ang2bohr = 1.8897261254578281

    # Add metal to the vdw_radii
    GridOptions = psiresp.grid.GridOptions()

    if metal_radius is not None and metal_name is not None:
        GridOptions.vdw_radii[metal_name] = metal_radius

    # create reorientations of the single conformer
    molecule_psiresp.generate_transformations(n_reorientations=n_reorientations)
    molecule_psiresp.generate_orientations(True)
    # GridOptions._generate_vdw_grid(np.array([atom.atomic_symbol for atom in site.atoms]), np.array([atom.coordinate.real for atom in site.atoms]))
    conf = molecule_psiresp.conformers[0]

    name = site.name

    # Iteration over the reorientations of the molecule
    for orient_idx in range(n_reorientations):

        # orient_idx=0
        orient = conf.orientations[orient_idx]
        # this does not work if already ortogonalized. it will break if file is reused.... :-(
        orient.compute_grid(GridOptions)
        print(orient.grid)

        # create grid for ORCA, which uses Bohr as length
        with open("grid.vpot.xyz", "w") as File:
            File.write(f"{len(orient.grid):}\n")
            for line in ang2bohr * orient.grid:
                File.write(f"{line[0]:} {line[1]:} {line[2]:}\n")

        ade.wrappers.ORCA.ORCA.execute = new_execute
        site.coordinates = orient.coordinates
        site.name = name + "_orient" + str(orient_idx)

        #method_keywords = ['B3LYP', '6-31G*', 'keepdens']
        #if metal_name is not None:
        #    if name_to_atomic_number[metal_name] > 30:


        site.single_point(method=method, keywords=keywords)

        # load esp from the file
        esp_name = site.name + '_sp_' + method.name + '.vpot.out'
        esp = np.loadtxt(esp_name, skiprows=1).T[3]

        # add esp and grid to the specific orientation
        molecule_psiresp.conformers[0].orientations[orient_idx].grid = orient.grid
        molecule_psiresp.conformers[0].orientations[orient_idx].esp = esp

    constraints = psiresp.ChargeConstraintOptions()
    if extra_atoms is not None:
        constrained_atoms = [index for index, _ in enumerate(site.atoms) if index in extra_atoms]

        logger.info(f"Constraining atoms {constrained_atoms:} Symbols: {[molecule_psiresp.atoms[a].symbol for a in constrained_atoms]:}")

        constraints.add_charge_sum_constraint_for_molecule(molecule_psiresp, charge=0, indices=constrained_atoms)
        logger.info(f"Constrains: {constraints.charge_sum_constraints:}")

    # run psiRESP to find the charges
    job = psiresp.Job(molecules=[molecule_psiresp], charge_constraints=constraints)
    normal_charges = job.compute_charges()

    resp_charges = normal_charges[0]
    logger.info(f"RESP before making extra atoms zero: {resp_charges:} Sum: {np.sum(resp_charges)}")

    # for some reason there is small charge left on residues
    # we calculatte mean of all, but not constrained atoms (we will set up them to zero)
    if extra_atoms is not None:
        mean_left_charge = (charge - np.sum(resp_charges) + np.sum(resp_charges[constrained_atoms])) / float(
            len(resp_charges) - len(constrained_atoms))
    else:
        mean_left_charge = (charge - np.sum(resp_charges)) / float(len(resp_charges))

    resp_charges += mean_left_charge

    if extra_atoms is not None:
        # we simply make the charges of the constraint 0, as the group is 0 (later on is easier to calculate difference)
        for idx in constrained_atoms:
            resp_charges[idx] = 0.0

    logger.info(f"RESP: {resp_charges:}, Sum: {np.sum(resp_charges):}")

    ade.wrappers.ORCA.ORCA.execute = old_execute
    os.chdir(here)
    return resp_charges

def calcualte_diffrence_and_symmetrize(metal_charge, unique_ligands_pattern, site_charges, ligand_charges):
    '''
    Calculates residual charges by subtracting charges of separate residues from the whole structure.
    Firstly the charges are symmetrized over the different ligands

    :param metal_charge: (int) charge of metal
    :param unique_ligands_pattern: (list(int)) ligand pattern, the same ligands get the same index
    :param site_charges: (list(float)) partial charges of the whole site
    :param ligand_charges: (list(float)) partial charges of separated ligands
    :return:
    '''
    logger.info(f"calcualte_diffrence() arguments: {metal_charge:}, {unique_ligands_pattern:}, {site_charges:}, {ligand_charges:}")
    # unconcatanatet to the split charges:
    n_lingads = [len(temp) for temp in ligand_charges]

    # number of atoms site
    site_atoms = [n_lingads[a] for a in unique_ligands_pattern]
    site_atoms = [1] + site_atoms
    split_site_charges = [np.array(site_charges[sum(site_atoms[:a]):sum(site_atoms[:a + 1])]) for a in
                          range(len(site_atoms))]
    logger.info(f"Charges of the whole site: {np.sum(np.concatenate(split_site_charges)):}")

    logger.info(f"Changing metal charge from {metal_charge} to {split_site_charges[0][0]}")
    split_site_charges[0][0] = split_site_charges[0][0] - metal_charge

    # we group the ligands and calculate mean value of charges
    for a in list(set(unique_ligands_pattern)):
        find = [c + 1 for c, b in enumerate(unique_ligands_pattern) if b == a]
        # +1 becasue of metal

        mean_split = np.mean(np.array([split_site_charges[idx] for idx in find]), axis=0)
        logger.info(f"Mean charge of unique ligands from whole site: {mean_split:}")
        for d in find:
            split_site_charges[d] = mean_split - ligand_charges[a]

    new_charges = np.concatenate(split_site_charges)
    logger.info(f"Check, sum of residual charge (should be ~0): {np.sum(new_charges):}")
    # Since this is residual charge, the sum of it should be below 0.1
    assert np.abs(np.sum(new_charges)) < 0.1

    return new_charges


def calculate_charges2(metal_name, metal_charge, filename, unique_ligands_pattern, ligand_charges_pattern, link_atoms,
                       additional_atoms, starting_index, metal_radius, mult=1):
    extra_atoms = link_atoms + additional_atoms

    unique_ligands_constraints = {}
    for key in list(set(unique_ligands_pattern)):
        unique_ligands_constraints[key] = []

    for extra_atom in extra_atoms:
        for idx, _ in enumerate(starting_index[:-1]):
            if starting_index[idx] < extra_atom < starting_index[idx + 1]:
                unique_ligands_constraints[unique_ligands_pattern[idx - 1]].append(
                    extra_atom - starting_index[idx])

    for idx, _ in enumerate(unique_ligands_constraints):
        unique_ligands_constraints[idx] = list(set(unique_ligands_constraints[idx]))

    ligand_charges = []
    for lig_idx in list(set(unique_ligands_pattern)):
        ligand_formal_charge = ligand_charges_pattern[list(unique_ligands_pattern).index(lig_idx)]
        charges = resp_orca(f"ligand_{lig_idx:d}.xyz", charge=ligand_formal_charge, opt=False, metal_name=None,
                            metal_radius=metal_radius, extra_atoms=unique_ligands_constraints[lig_idx])
        ligand_charges.append(charges)

    site_charges = resp_orca('saturated_template.xyz', charge=metal_charge + np.sum(ligand_charges_pattern), opt=False,
                             metal_name=metal_name, metal_radius=metal_radius, mult=mult, extra_atoms=extra_atoms)

    print("ligand charges", len(ligand_charges), ligand_charges)
    print("ligand        ", len(site_charges), site_charges)

    new_site_charges = calcualte_diffrence_and_symmetrize(metal_charge, unique_ligands_pattern, site_charges,
                                                          ligand_charges)

    return new_site_charges

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-charge", help="Charge")
    parser.add_argument("-metal_name", help="Name of the metal")
    parser.add_argument("-vdw_data_name", default='uff', help="vdw_data_name")
    parser.add_argument("-mult", default=1, help="multiplicity")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    # if args.o is not None:
    charges = resp_orca(args.f, charge=args.charge, metal_name=args.metal_name, metal_radius=args.metal_radius, mult=args.mult, extra_atoms=args.extra_atoms)