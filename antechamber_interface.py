import parmed as pmd

from tempfile import mkdtemp
import os
import shutil
from subprocess import Popen, DEVNULL
from log import logger




def antechamber(pdbfile, charge, output, verbose=False):

    def run_external(command, assertion=None):
        with open("output.txt", 'w') as output_file:
            process = Popen(command.split(), stdout=output_file, stderr=DEVNULL)
            process.wait()
        logger.info("COMMAND")
        logger.info(command)

        if assertion is not None:
            with open("output.txt") as File:
                text = File.read()
                logger.info(text)
        #        assert assertion in text


    pwd = os.getcwd()
    tmpdir_path = mkdtemp()
    logger.info(f"[ ] Current path: {pwd:s}") if verbose else None
    logger.info(f"[ ] Going to temporary path {tmpdir_path:s}") if verbose else None
    os.chdir(tmpdir_path)
    shutil.copyfile(f'{pwd:s}/{pdbfile:s}', 'temp.pdb')

    logger.info("[ ] Interfacing antechamber") if verbose else None
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

    logger.info("    - Charge of the molecule is", charge)  if verbose else None
    run_external("antechamber -i temp.pdb -fi pdb -o temp.mol2 -fo mol2 -c bcc -s 2 -nc "+str(charge), assertion="Errors = 0")
    run_external("parmchk2 -i temp.mol2 -f mol2 -o temp.frcmod")
    run_external("tleap -f tleap.in")

    parm = pmd.load_file('temp.prmtop', 'temp.inpcrd')
    parm.save('topol.top', format='gromacs')
    shutil.copyfile('topol.top', f'{pwd:s}/{output:s}')
    logger.info("[ ] Molecule parametrized")  if verbose else None

    os.chdir(pwd)
    shutil.rmtree(tmpdir_path)

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="pdb file")
    parser.add_argument("-charge", default=0, help="Topology of the linker ") # TODO
    parser.add_argument("-o", default='topol.top', help="Output coordination file")
    parser.add_argument("-v", action='store_true', dest='v', help="Be loud")
    parser.set_defaults(v=False)
    return parser.parse_args()


# MAKE it less louder
if __name__ == '__main__':
    args = get_args()
    antechamber(args.f, args.charge, args.o, args.v)