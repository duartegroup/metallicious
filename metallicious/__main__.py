from metallicious.parametrize_new_sites import supramolecular_structure
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-p", help="Metaloorganic topology (*top, *prmtop, etc. all supported by ParmEd)")

    #parser.add_argument("-o", help="Clean structure")

    parser.add_argument("-metal_and_charges",nargs='+', help="Metal names and charges (optionally, multiplicity, when parametrization needed), example: Pd 2 1 Ru 2 1")
    parser.add_argument("-keywords", help="keywords for QM", nargs='+')
    parser.add_argument("-vdw_type", default='merz-opc',
                        help="Type of parameters for VdW (available: uff, merz-tip3p, merz-opc3, merz-spc/e, merz-tip3p-fb, merz-opc, merz-tip4p-fb, merz-tip4-ew, zhang-tip3p, zhang-opc3, zhang-spc/e, zhang-spc/eb, zhang-tip3p-fb, zhang-opc, zhang-tip4p/2005, zhang-tip4p-d, zhang-tip4p-fb, zhang-tip4p-ew")
    parser.add_argument("-truncate", help="Truncatation scheme (available: None, dihedral, angle, bond)", default=None)
    parser.add_argument("-improper_metal", action='store_true', default=False,
                        help="Calculate improper dihedral of the metal-aromatic (default:False)")
    parser.add_argument("-donors", nargs='+', default=['N', 'S', 'O'],
                        help="Donors from the connected ligands, usually electronegative atom, such as N, S, O, but sometimes metal is connected to carbon", )
    return parser.parse_args()

def main():
    args = get_args()

    if args.f is None:
        print("The coordination file is required")
        raise
    else:
        filename = args.f

    if args.f is None:
        print("The topology file is required")
        print("Simple GAFF parametrization can be done by XXX") # TODO
        raise
    else:
        topol = args.p

    is_metal_charges_mult = False
    if len(args.metal_and_charges) % 2==0 or len(args.metal_and_charges) % 3==0:
        if len(args.metal_and_charges)>2:
            # if args.metal_and_charges[2] is a number: # TODO
            is_metal_charges_mult = True
            None
    else:
        print("Not correct format of metal_and_charges")
        return False

    metal_charge_mult = None
    metal_charges = None

    if is_metal_charges_mult:
        if len(args.metal_and_charges) % 3==0:
            metal_charge_mult = {}
            for idx in range(int(len(args.metal_and_charges)/3)):
                metal_charge_mult[args.metal_and_charges[idx*3]] = (int(args.metal_and_charges[idx*3+1]), int(args.metal_and_charges[idx*3+2]))
        else:
            print("Not correct format of metal_charges")
            raise
    else:
        if len(args.metal_and_charges) % 2==0:
            metal_charges = {}
            for idx in range(int(len(args.metal_and_charges)/2)):
                metal_charges[args.metal_and_charges[idx*2]] = int(args.metal_and_charges[idx*2+1])

    if args.keywords is None:
        keywords = ['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq']
    else:
        keywords = args.keywords

    library_path = f'{os.path.dirname(__file__):s}/library/'

    vdw_type = args.vdw_type

    improper_metal = args.improper_metal
    donors = args.donors
    truncation_scheme = args.truncate

    cage = supramolecular_structure(filename=filename, metal_charge_mult=metal_charge_mult,
                                    metal_charges=metal_charges, vdw_type=vdw_type, topol=topol, keywords=keywords,
                                    improper_metal=improper_metal, donors=donors, library_path=library_path,
                                    truncation_scheme=truncation_scheme)
    cage.parametrize()

if __name__ == '__main__':
    main()
