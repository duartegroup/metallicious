import argparse


multiwfn_path = "/home/fd05/fd/chem1540/Programs/Multiwfn/Multiwfn"

import autode as ade
from autode.wrappers.ORCA import *
from subprocess import Popen, PIPE, STDOUT
import MDAnalysis
import psiresp
import qcelemental as qcel

import os

try:
    from cgbind2pmd.data import vdw_data
    from cgbind2pmd.data import name_to_atomic_number
except:
    from data import vdw_data
    from data import name_to_atomic_number

def perform_resp_with_multiwfn(molden_input):
    '''
    Depracticed function as Multiwfn is hard to use on cluster and it cannot be installed from conda

    :param molden_input:
    :return:
    '''
    File = open("molden_input.txt", "w")
    File.write("7\n18\n1\n\ny\nq\n")
    File.close()

    '''
    with open("molden_input.txt") as infile:
        process = Popen([multiwfn_path, molden_input, "-silent"], stdin=infile)
        process.wait()
    '''
    filename_core = molden_input.replace(".input", "")

    File = open(filename_core + ".chg")
    text = File.read()
    File.close()
    charges = list(map(float, text.split()[4::5]))
    return charges

def perform_resp_with_psiresp(molden_input):
    None


def old_execute(self, calc):
    @work_in_tmp_dir(filenames_to_copy=calc.input.filenames,
                     kept_file_exts=('.out', '.hess', '.xyz', '.inp', '.pc'))
    def execute_orca():
        # Run the calculations
        run_external(params=[calc.method.path, calc.input.filename],
                     output_filename=calc.output.filename)

    execute_orca()
    return None





#@check_sufficient_memory
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
        process = Popen(params, stdout=output_file, stderr=PIPE,  close_fds=True)

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
    @work_in_tmp_dir(filenames_to_copy=calc.input.filenames+["grid.vpot.xyz"], kept_file_exts=('.out', '.hess', '.xyz', '.inp', '.pc', '.input', '.densities', '.gbw'))
    def execute_orca():
        # print([calc.method.path, calc.input.filename], calc.output.filename)
        # Run the calculations
        run_external(params=[calc.method.path, calc.input.filename],
                     output_filename=calc.output.filename)
        # Created molden file
        #filename_core = calc.input.filename.replace(".inp", "")
        #run_external(params=[calc.method.path + "_2mkl", filename_core, "-molden"],
        #                     output_filename=filename_core + ".input")
        print(os.listdir())
        #orca_vpot filename.gbw filename.scfp filename.vpot.xyz filename.vpot.out
        filename_core = calc.input.filename.replace(".inp", "")
        run_external(params=[calc.method.path+"_vpot", filename_core+".gbw", filename_core+".scfp", "grid.vpot.xyz", filename_core+".vpot.out"],
                     output_filename=filename_core+"_vpot.out")
    execute_orca()
    return None

def resp_orca(filename, charge=0, opt=True, metal_name=None, vdw_data_name='uff', n_reorientations = 1, mult=1, extra_atoms=None):
    '''

    import sys
    sys.path.insert(0, '/u/fd/chem1540/github/')
    from cgbind2pmd.charges import *

    filename='site.xyz'
    charge=2
    vdw_data_name='uff'
    n_reorientations = 1
    metal_name='Fe'

    '''
    # TODO you cannot restart resp! you need to delete folder

    method = ade.methods.ORCA()
    site = ade.Molecule(filename, charge=charge, mult=mult)
    print(filename)
    print(site, site.atoms)

    molecule_qcel = qcel.models.Molecule.from_file(filename, molecular_charge=charge)
    molecule_psiresp = psiresp.Molecule(qcmol=molecule_qcel, charge=charge)
    print("molecule_qcel", molecule_qcel)
    print("molecule_psiresp", molecule_psiresp)

    here = os.getcwd()
    os.system("mkdir resp")
    os.chdir("resp")

    if opt:
        ade.wrappers.ORCA.ORCA.execute = old_execute
        site.optimise(method=method)


    data = vdw_data[vdw_data_name]
    print("data", data)

    ang2bohr = 1.8897261254578281

    # Add metal to the vdw_radii
    GridOptions = psiresp.grid.GridOptions()

    if metal_name is not None:
        GridOptions.vdw_radii[metal_name] = data[metal_name][1] # psiresp uses R distance rather then R/2 TODO check for sure

    # create reorientations of the single conformer
    molecule_psiresp.generate_transformations(n_reorientations=n_reorientations)
    molecule_psiresp.generate_orientations()
    #GridOptions._generate_vdw_grid(np.array([atom.atomic_symbol for atom in site.atoms]), np.array([atom.coordinate.real for atom in site.atoms]))
    conf = molecule_psiresp.conformers[0]

    name = site.name

    # Iteretation over the reorientations of the molecule
    for orient_idx in range(n_reorientations):

        #orient_idx=0
        orient = conf.orientations[orient_idx]
        orient.compute_grid(GridOptions) # this does not work if already ortognaized... it will break if at exactly 0,0,0 TODO
        print(orient.grid)

        # create frid for ORCA, which uses Bohr as lenght
        with open("grid.vpot.xyz", "w") as File:
            File.write(f"{len(orient.grid):}\n")
            for line in ang2bohr * orient.grid:
                File.write(f"{line[0]:} {line[1]:} {line[2]:}\n")

        ade.wrappers.ORCA.ORCA.execute = new_execute
        site.coordinates = orient.coordinates
        site.name = name + "_orient" + str(orient_idx)

        #basis_on_the_metal = f'\n%basis\nnewgto {metal_name:s} "LANL2DZ" end\nend\n'

        method_keywords = ['B3LYP', '6-31G*', 'keepdens']
        if metal_name is not None:
            if name_to_atomic_number[metal_name] > 30:
                method_keywords = ['PBE0', 'def2-SVP', 'keepdens']
        site.single_point(method=method, keywords=method_keywords)


        # load esp from the file
        esp_name = site.name + '_sp_' + method.name + '.vpot.out'
        esp = np.loadtxt(esp_name, skiprows=1).T[3]

        # add esp and grid to the specific orientation
        molecule_psiresp.conformers[0].orientations[orient_idx].grid = orient.grid
        molecule_psiresp.conformers[0].orientations[orient_idx].esp = esp

    constraints = psiresp.ChargeConstraintOptions()
    if extra_atoms is not None:
        constrained_atoms = [index for index, _ in enumerate(site.atoms) if index in extra_atoms]

        print("Constraining atoms", constrained_atoms, "Symbols:", [molecule_psiresp.atoms[a].symbol for a in constrained_atoms])
        print("Constrains:", constraints.charge_sum_constraints)

        constraints.add_charge_sum_constraint_for_molecule(molecule_psiresp, charge=0, indices=constrained_atoms)
        print("Constrains:", constraints.charge_sum_constraints)

    # run psiRESP to find the charges
    job = psiresp.Job(molecules=[molecule_psiresp],  charge_constraints=constraints)
    normal_charges = job.compute_charges()

    print("RESP charges", normal_charges[0])
    if extra_atoms is not None:
        print("Constrained charges: ", normal_charges[0][constrained_atoms], "Sum:", np.sum(normal_charges[0][constrained_atoms]))

    resp_charges = normal_charges[0]
    print("before RESP:", resp_charges, "Sum:", np.sum(resp_charges))

    #for some reason there is small charge left on this residue
    # we calculatte mean of all, but not constrained atoms (we will set up them to zero)
    if extra_atoms is not None:
        mean_left_charge = (charge-np.sum(resp_charges)+np.sum(resp_charges[constrained_atoms]))/float(len(resp_charges)-len(constrained_atoms))
    else:
        mean_left_charge = (charge - np.sum(resp_charges))/float(len(resp_charges))
    print("mean_left_charge", mean_left_charge)
    resp_charges += mean_left_charge

    if extra_atoms is not None:
        # we simply make the charges of the constrain 0, as the group is 0 (later on is easier to calculate diffrence)
        for idx in constrained_atoms:
            resp_charges[idx] = 0.0

    print("RESP:", resp_charges, "Sum:", np.sum(resp_charges))


    os.chdir(here)
    return resp_charges


def resp_psi4(filename, charge=0, opt=True, metal_name=None):
    try:
        import psi4
        import resp
    except:
        print("you need psi4!")
        return False

    with open(filename) as f:
        xyz = f.read()

    mol = psi4.geometry(xyz)
    mol.update_geometry()
    mol.set_molecular_charge(charge)

    basis = '6-31G*'
    if metal_name is not None:
        if name_to_atomic_number[metal_name]>36:
            basis='lanl2dz'

    if opt:
        psi4.optimize(f'pbe0/def2-svp')

    options = {'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0],
               'VDW_POINT_DENSITY': 1.0,
               'RESP_A': 0.0005,
               'RESP_B': 0.1,
               'METHOD_ESP': 'B3LYP',
               'BASIS_ESP': basis}

    # Call for first stage fit
    charges1 = resp.resp([mol], options)
    print('Electrostatic Potential Charges')
    print(charges1[0])
    print('Restrained Electrostatic Potential Charges')
    print(charges1[1])

    # Change the value of the RESP parameter A
    options['RESP_A'] = 0.001

    # Add constraint for atoms fixed in second stage fit
    options['grid'] = ['1_%s_grid.dat' % mol.name()]
    options['esp'] = ['1_%s_grid_esp.dat' % mol.name()]

    # Call for second stage fit
    charges2 = resp.resp([mol], options)

    # Get RESP charges
    print("\nStage Two:\n")
    print('RESP Charges')
    print(charges2[1])

    return charges2[1]


def calcualte_diffrence(metal_charge, unique_ligands_pattern, site_charges, ligand_charges):
    print("calcualte_diffrence", metal_charge, unique_ligands_pattern, site_charges, ligand_charges)
    #unconcatanatet to the split charges:
    # atoms in ligands
    n_lingads = [len(temp) for temp in ligand_charges]

    # number of atoms site
    site_atoms = [n_lingads[a] for a in unique_ligands_pattern]
    site_atoms = [1] + site_atoms
    split_site_charges = [np.array(site_charges[sum(site_atoms[:a]):sum(site_atoms[:a+1])]) for a in range(len(site_atoms)) ]

    print(split_site_charges[0][0])
    print(metal_charge)

    split_site_charges[0][0] = split_site_charges[0][0] - metal_charge

    # we group the ligands and calculate mean value of charges
    for a in list(set(unique_ligands_pattern)):
        find = [c + 1 for c, b in enumerate(unique_ligands_pattern) if b == a]
        # +1 becasue of metal
        print("Grouped ligands", find)
        mean_split = np.mean(np.array(split_site_charges)[find])
        print("mean_split", mean_split)
        for d in find:
            split_site_charges[d] = mean_split - ligand_charges[a]

    new_charges = np.concatenate(split_site_charges)
    return new_charges


def calculate_charges(metal_charge, metal_name, vdw_data_name, mult=1):

    File = open("INFO.dat")
    text = File.read()
    File.close()

    n_sites = text.count("ligand_pattern")
    print("Number of sites", n_sites)

    unique_ligands_pattern = None
    all_sites_charges = []

    n_site = 0
    for line in text.splitlines():
        if "ligand_pattern:" in line:
            unique_ligands_pattern = list(map(int, line[15:].split(',')))

        if "link_atoms:" in line:
            if ',' in line:
                link_atoms = list(map(int, line[11:].split(',')))
            else:
                link_atoms = []

        if "extra_atoms:" in line:
            if ',' in line:
                additional_atoms = list(map(int, line[12:].split(',')))
            else:
                additional_atoms = []

        if "starting_index:" in line:
            starting_index = list(map(int, line[15:].split(',')))

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
                charges = resp_orca(f"ligand{n_site:d}_{lig_idx:d}.xyz", opt=False, metal_name=metal_name, vdw_data_name=vdw_data_name, extra_atoms=unique_ligands_constraints[lig_idx])
                ligand_charges.append(charges)

            site_charges = resp_orca(f"site{n_site:d}.xyz", charge=metal_charge, opt=False, metal_name=metal_name, vdw_data_name=vdw_data_name, mult=mult, extra_atoms=extra_atoms)

            print("ligand charges", len(ligand_charges), ligand_charges)
            print("ligand        ", len(site_charges), site_charges)

            new_site_charges = calcualte_diffrence(metal_charge, unique_ligands_pattern, site_charges, ligand_charges)
            all_sites_charges.append(new_site_charges)
            n_site+=1

    print("all_sites", all_sites_charges)

    return all_sites_charges


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
    #if args.o is not None:
    print(calculate_charges(float(args.charge), args.metal_name, args.vdw_data_name), mult=args.mult)

