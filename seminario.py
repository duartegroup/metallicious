path_to_mod_seminario = "/u/fd/chem1540/github/ModSeminario_Py/Python_Modified_Seminario_Method"


import autode as ade
method = ade.methods.XTB()
import os
import re
from subprocess import Popen, DEVNULL


name_to_atomic_number = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                         'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
                         'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
                         'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38,
                         'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47,
                         'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
                         'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
                         'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
                         'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83,
                         'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
                         'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101,
                         'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
                         'Ds': 110, 'Rg': 111, 'Uub': 112}


def orca_to_fchk(filename="site_opt_orca.hess", output_fchk="lig.fchk"):
    File = open(filename)
    text = File.read()
    File.close()

    temp = text[text.find('$atoms'):text.find("$actual_temperature")]
    N_atoms = int(temp.splitlines()[1])

    gaussian_coord_string = ''
    gaussian_atomic_number_string = ''

    c = 0
    d = 0
    for line in temp.splitlines()[2:]:
        print(line)

        if len(line.split()) > 0:
            for a in [2, 3, 4]:

                gaussian_coord_string += f" {float(line.split()[a]): .8e}"
                c += 1
                if c % 5 == 0:
                    gaussian_coord_string += '\n'
                    c = 0
            gaussian_atomic_number_string += f" {name_to_atomic_number[line.split()[0].title()]: 11d}"
            d += 1

            if d % 5 == 0:
                gaussian_atomic_number_string += '\n'
                d = 0

    temp = text[text.find('$hessian'):text.find("$vibrational_frequencies")]
    N = int(temp.splitlines()[1])

    hess_string = ['' for a in range(N)]

    for a in range(int(N / 5) + 1):
        temp2 = temp.splitlines()[3 + (1 + N) * (a):2 + (1 + N) * (a + 1)]
        for b in range(N):
            hess_string[b] += temp2[b][6:]

    gaussian_hess_string = ''

    c = 0
    for a in range(N):
        for b in range(a + 1):
            gaussian_hess_string += f" {float(hess_string[b].split()[a]): .8e}"
            c += 1
            if c % 5 == 0:
                gaussian_hess_string += '\n'
                c = 0

    File = open(output_fchk, "w")
    print(f"Atomic numbers                             I   N={N_atoms:12d}", file=File)
    print(gaussian_atomic_number_string, file=File)
    print("Nuclear charges        \n", file=File)
    print(f"Current cartesian coordinates              R   N={N:12d}", file=File)
    print(gaussian_coord_string, file=File)
    print("Force Field            \n", file=File)
    print(f"Cartesian Force Constants                  R   N={N:12d}", file=File)
    print(gaussian_hess_string, file=File)
    print("Dipole Moment          \n", file=File)
    File.close()


def orca_to_log(filename="site_opt_orca.out", log_output="lig.log"):
    File = open(filename)
    text = File.read()
    File.close()

    start_line = text.find("        Definition                    Value    dE/dq     Step     New-Value")
    n_bonds = 1
    n_angles = 1

    File = open(log_output, "w")

    print(" --------------------------                            --------------------------", file=File)
    print(" ! Name  Definition              Value          Derivative Info.                !", file=File)
    print(" --------------------------------------------------------------------------------", file=File)

    for line in text[start_line:].splitlines()[2:]:
        if line == "    ----------------------------------------------------------------------------":
            break

        if line[8] == 'B':

            match = re.search("\([a-zA-Z]+\s*(\d+),[a-zA-Z]+\s*(\d+)\)", line)
            idx1 = int(match.group(1)) + 1
            idx2 = int(match.group(2)) + 1

            number_of_spaces = 22 - len(f"R({idx1:d},{idx2:d})")
            value = float(line[line.find(")"):].split()[1])

            print(
                f" ! R{n_bonds:<4d} R({idx1:d},{idx2:d}){number_of_spaces * ' ':s}{0.0: 8.4f}         -DE/DX =    0.0                 !",
                file=File)
            n_bonds += 1

        elif line[8] == 'A':
            match = re.search("A\([a-zA-Z]+\s*(\d+),[a-zA-Z]+\s*(\d+),[a-zA-Z]+\s*(\d+)\)", line)
            idx1 = int(match.group(1)) + 1
            idx2 = int(match.group(2)) + 1
            idx3 = int(match.group(3)) + 1

            number_of_spaces = 22 - len(f"A({idx1:d},{idx2:d},{idx3:d})")

            value = float(line[line.find(")"):].split()[1])

            print(
                f" ! A{n_angles:<4d} A({idx1:d},{idx2:d},{idx3:d}){number_of_spaces * ' ':s}{value:8.4f}         -DE/DX =    0.0                 !",
                file=File)

            n_angles += 1
    print(" ! D1    D(1,2,3,4)              0.0000         -DE/DX =    0.0                 !", file=File)
    print(" --------------------------------------------------------------------------------\n", file=File)
    File.close()


def frequencies(filename, charge = 0):
    site = ade.Molecule(filename)
    site.charge = charge
    site.optimise(method=method, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'])

def seminario(filename, metal_charge):
    here = os.getcwd()
    frequencies(filename, metal_charge)
    orca_to_fchk("site_opt_orca.hess")
    orca_to_log("site_opt_orca.out")

    os.chdir(path_to_mod_seminario)
    # command = f'python modified_Seminario_method.py {here:s}/ {here:s}/ 0.957'
    command = f'python modified_Seminario_method.py {here:s}/ {here:s}/ 1.000'
    process = Popen(command.split(), stdout=DEVNULL, stderr=DEVNULL)
    process.wait()
    os.chdir(here)





