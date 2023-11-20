import parmed as pmd

from tempfile import mkdtemp
import os
import shutil

from subprocess import Popen, DEVNULL
import MDAnalysis

# try:
from metallicious.log import logger
# except:
#     from log import logger
#     from utils import mdanalysis_to_rdkit



import rdkit
import numpy as np

def check_antechamber_if_available():
    '''
    Checks if antechamber (ambertools) is installed)
    :return:
    '''
    return shutil.which('antechamber') is not None


def antechamber(pdbfile,output, charge=None):
    '''
    Simple antechamber interface.  It runs with GAFF2. Sometimes it gives charge which is not equal to formal charge:
    we simply subtract the diffrence over all atoms

    :param pdbfile: (string) pdb file
    :param output:  (string) top output file
    :param charge: (int) formal charge of a molecule
    :return: None
    '''

    def run_external(command, assertion=None):
        with open("output.txt", 'w') as output_file:
            process = Popen(command.split(), stdout=output_file, stderr=DEVNULL)
            process.wait()
        logger.info(command)

        if assertion is not None:
            with open("output.txt") as File:
                text = File.read()
                logger.info(text)
        #        assert assertion in text


    pwd = os.getcwd()
    tmpdir_path = mkdtemp()
    logger.info(f"[ ] Current path: {pwd:s}")
    logger.info(f"[ ] Going to temporary path {tmpdir_path:s}")

    syst = MDAnalysis.Universe(pdbfile)
    os.chdir(tmpdir_path)

    if len(syst.atoms.names)!=len(set(syst.atoms.names)):
        for a, atom in enumerate(syst.atoms):
            atom.name = atom.name +str(a)
    syst.atoms.write('temp.pdb')

    if charge is None:
        #mol = mdanalysis_to_rdkit(syst.atoms)
        mol = syst.atoms.convert_to("RDKIT")
        charge = rdkit.Chem.GetFormalCharge(mol)
        logger.info(f"    Guessing charge: {charge:}")

    logger.info("[ ] Interfacing antechamber")
    File = open("tleap.in", "w")
    File.write('''
    source leaprc.gaff\n
    TEMP = loadmol2 temp.mol2\n
    check TEMP \n
    loadamberparams temp.frcmod\n
    saveoff TEMP sus.lib \n
    saveamberparm TEMP temp.prmtop temp.inpcrd\n
    quit\n''')
    File.close()

    logger.info(f"    - Charge of the molecule is {charge:}")
    #sqm_key = "grms_tol=0.005,scfconv=1.d-8,ndiis_attempts=700,maxcyc=0,"
    #run_external(f'antechamber -i temp.pdb -fi pdb -o temp.mol2 -fo mol2 -c bcc -at gaff2 -s 2 -nc {charge:} -ek "{sqm_key:}"', assertion="Errors = 0")
    run_external(
        f'antechamber -i temp.pdb -fi pdb -o temp.mol2 -fo mol2 -c bcc -at gaff2 -s 2 -nc {charge:}',
        assertion="Errors = 0")
    run_external("parmchk2 -i temp.mol2 -f mol2 -o temp.frcmod")
    run_external("tleap -f tleap.in")

    parm = pmd.load_file('temp.prmtop', 'temp.inpcrd')
    parm.save('topol.top', format='gromacs')


    neutralize_charge('topol.top', 'topol2.top', charge=charge)

    shutil.copyfile('topol2.top', f'{pwd:s}/{output:s}')
    logger.info("[ ] Molecule parametrized")

    os.chdir(pwd)
    shutil.rmtree(tmpdir_path)

    syst.atoms.write(f'{pwd:s}/{output[:-4]:s}.pdb')


import argparse


def neutralize_charge(file_name, output, charge=0):
    '''
    Normalizes charge present in *top file to specified charge.
    It is to prevent antechamber from small diffrences of charge from expected one.
    :param file_name: (string) input topology file in gromacs format (*top)
    :param output: (string) output topology file in gromacs format (*top)
    :param charge: (int) charge
    :return: None
    '''
    File = open(file_name, "r")
    text = File.read()
    File.close()
    moleculetype = ""

    if text.count("[ moleculetype ]") > 1:
        logger.info("ERROR: More than one moleculetype")
    elif ("[ system ]" in text):
        moleculetype = text[text.find(r"[ moleculetype ]"):text.find("[ system ]")]
    else:
        moleculetype = text[text.find("[ moleculetype ]"):]

    top = output
    File = open(top, "w")
    File.write(text[:text.find("[ moleculetype ]")])

    atom = True
    in_atoms = False

    # Correction for charge, should be inside the parmetrization
    sum_of_charge = - float(charge)  # we start with the formal charge
    number_of_atoms = 0

    for line in moleculetype.splitlines():
        if in_atoms and len(line) > 0 and line[0] != ";":
            if len(line.split()) > 7:
                sum_of_charge += float(line.split()[6])
                number_of_atoms += 1

        if "[ atoms ]" in line:
            in_atoms = True
        if "[ bonds ]" in line:
            in_atoms = False

    if abs(sum_of_charge) > 0.00001:
        fract = sum_of_charge / number_of_atoms
        fract_rest = sum_of_charge - number_of_atoms * np.round(fract, 6)
        every_nth = int(number_of_atoms / (np.abs(fract_rest) / 0.000001))
        qtot = 0.0
        logger.info(f"   [ ] Rounding up charges to make molecule neutral, every atom gets additional {-fract:}")
        logger.info(f"       ,every atom will be assigned 0.000001 each {fract_rest:}")
        atom_nr = 0
        for line in moleculetype.splitlines():
            # print(np.round(fract + 0.000001*np.random.random(),6)) #
            if in_atoms and len(line) > 0 and line[0] != ";" and len(line.split()) > 6:
                line = line[:line.find(';')]
                charge = float(line.split()[6]) - np.round(fract + np.sign(fract_rest) * 0.000001 * (
                        atom_nr % every_nth == 0 and every_nth * int(
                    (np.abs(fract_rest) / 0.000001)) > atom_nr), 6)
                qtot += charge  # w

                print("{:>6s}{:>11s}{:>7s}{:>8s}{:>6s}{:>7s}{:>10f}{:>11s} ; qtot {: 6f}".format(*line.split()[:6],
                                                                                                 charge,
                                                                                                 line.split()[7],
                                                                                                 qtot),
                      file=File)  # ,[str(np.round(float(line.split()[6])-fract,6))],line.split()[7]))                else:
                atom_nr += 1
            else:
                print(line, file=File)
                # print(qtot)

            if "[ atoms ]" in line:
                in_atoms = True
            if "[ bonds ]" in line:
                in_atoms = False
    else:
        for line in moleculetype.splitlines():
            print(line, file=File)

    File.close()


def smiles_to_mol(smiles, outfile='rdkit.pdb'):
    '''
    Converts SMILES string into pdb file
    :param smiles: (string) SMILES string
    :param outfile: (string) output PDB file
    :return:
    '''

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    embeded = AllChem.EmbedMolecule(mol)
    if embeded == -1:  # if does not work use just random coords
        AllChem.EmbedMolecule(mol, useRandomCoords=True)

    Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, outfile)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="pdb file")
    parser.add_argument("-smi", help="smiles")
    parser.add_argument("-charge", default=0, help="Topology of the linker ")
    parser.add_argument("-o", default='topol.top', help="Output coordination file")
    parser.add_argument("-v", action='store_true', dest='v', help="Be loud")
    parser.set_defaults(v=False)
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    if args.f is not None:
        antechamber(args.f,  args.o,args.charge, args.v)
    elif args.smi is not None:
        smiles_to_mol(args.smi, 'rdkit.pdb')
        antechamber('rdkit.pdb', args.o, args.charge, args.v)
    else:
        print("Provide input file")