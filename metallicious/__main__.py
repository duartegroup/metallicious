from metallicious.parametrize_new_sites import supramolecular_structure
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structure (*.gro, *.pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-p", help="Metaloorganic topology of the whole structure (*.top, *.prmtop, etc. all supported by ParmEd)", default=False)

    parser.add_argument("-of", help="Output metaloorganic structure (*.gro, *.pdb, etc. supported by MDAnalysis)", default='out.pdb')
    parser.add_argument("-op", help="Output metaloorganic topology (*.top, *.prmtop, etc. supported by ParmEd)", default='out.top')

    parser.add_argument("-metal_and_charges",nargs='+', help="Metal names and charges (optionally, multiplicity, when parametrization needed), for example: Pd 2 1 Ru 2 1")
    parser.add_argument("-keywords", help="autodE keywords for QM calculations (default: PBE0 D3BJ def2-SVP tightOPT freq)", nargs='+')
    parser.add_argument("-vdw_type", default='merz-opc',
                        help="Type of parameters for Lennard-Jones paramters (default: merz-opc; available: uff, merz-tip3p, merz-opc3, merz-spc/e, merz-tip3p-fb, merz-opc, merz-tip4p-fb, merz-tip4-ew, zhang-tip3p, zhang-opc3, zhang-spc/e, zhang-spc/eb, zhang-tip3p-fb, zhang-opc, zhang-tip4p/2005, zhang-tip4p-d, zhang-tip4p-fb, zhang-tip4p-ew")
    parser.add_argument("-truncate", help="Truncation scheme (defualt: None; available: None, dihedral, angle, bond)", default=None)
    parser.add_argument("-improper_metal", action='store_true', default=False,
                        help="Calculate improper dihedral of the metal-aromatic (default: False)")
    parser.add_argument("-donors", nargs='+', default=['N', 'S', 'O'],
                        help = "Donors from the connected ligands, usually electronegative atom, such as N, S, O, but sometimes metal is connected to carbon (default: N S O)")
    parser.add_argument("-prepare_topol", action='store_true', default=False, help="Prepare initial topology using GAFF")
    parser.add_argument("-linker_topol", default=None, help="Linker force-field (topology) parameters")
    return parser.parse_args()

def main():
    args = get_args()

    if args.f is None:
        print("The coordination file is required, see: metallicious -h")
        return False
    else:
        filename = args.f

    if args.f is None:
        print("The topology file is required, see metallicious -h")
        print("Simple GAFF parametrization can be done by -prepare_topol")
        return False
    else:
        topol = args.p

    is_metal_charges_mult = False
    if len(args.metal_and_charges) % 2==0 or len(args.metal_and_charges) % 3==0:
        if len(args.metal_and_charges)>2:
            if args.metal_and_charges[2].isnumeric():
                is_metal_charges_mult = True
    else:
        print("Not correct format of metal_and_charges, see: metallicious -h")
        return False

    metal_charge_mult = None
    metal_charges = None

    if is_metal_charges_mult:
        if len(args.metal_and_charges) % 3==0:
            metal_charge_mult = {}
            for idx in range(int(len(args.metal_and_charges)/3)):
                metal_charge_mult[args.metal_and_charges[idx*3]] = (int(args.metal_and_charges[idx*3+1]), int(args.metal_and_charges[idx*3+2]))
        else:
            print("Not correct format of metal_charges, see: metallicious -h")
            return False
    else:
        if len(args.metal_and_charges) % 2==0:
            metal_charges = {}
            for idx in range(int(len(args.metal_and_charges)/2)):
                metal_charges[args.metal_and_charges[idx*2]] = int(args.metal_and_charges[idx*2+1])
        else:
            print("Not correct format of metal_charges, see: metallicious -h")
            return False


    if args.keywords is None:
        keywords = ['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq']
    else:
        keywords = args.keywords


    vdw_type = args.vdw_type

    improper_metal = args.improper_metal
    donors = args.donors
    truncation_scheme = args.truncate

    cage = supramolecular_structure(filename=filename, metal_charge_mult=metal_charge_mult,
                                    metal_charges=metal_charges, vdw_type=vdw_type, topol=topol, keywords=keywords,
                                    improper_metal=improper_metal, donors=donors,
                                    truncation_scheme=truncation_scheme)

    if args.prepare_topol:
        cage.prepare_initial_topology(homoleptic_ligand_topol=args.linker_topol)

    cage.parametrize(args.of, args.op)



if __name__ == '__main__':
    main()
