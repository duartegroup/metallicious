import parmed as pmd

from tempfile import mkdtemp
import os
import shutil

def antechamber(pdbfile, charge, output):
    print("[ ] Interfacing antechamber")
    pwd = os.getcwd()
    tmpdir_path = mkdtemp()
    os.chdir(tmpdir_path)
    shutil.copyfile(f'{pwd:s}/{pdbfile:s}', 'temp.pdb')

    File = open("tleap.in", "w")
    File.write("source leaprc.gaff\n")
    File.write("TEMP = loadmol2 temp.mol2 \n")
    File.write("check TEMP \n")
    File.write("loadamberparams temp.frcmod\n")
    File.write("saveoff TEMP sus.lib \n")
    File.write("saveamberparm TEMP temp.prmtop temp.inpcrd\n")
    File.write("quit\n")
    File.close()

    print("    - Charge of the molecule is", charge)
    os.system("antechamber -i temp.pdb -fi pdb -o temp.mol2 -fo mol2 -c bcc -s 2 -nc "+str(charge))
    os.system("parmchk2 -i temp.mol2 -f mol2 -o temp.frcmod")
    os.system("tleap -f tleap.in")

    parm = pmd.load_file('temp.prmtop', 'temp.inpcrd')
    parm.save('topol.top', format='gromacs')
    shutil.copyfile('topol.top', f'{pwd:s}/{output:s}')
    print("[ ] Molecule parametrized")
    os.chdir(pwd)
    shutil.rmtree(tmpdir_path)
