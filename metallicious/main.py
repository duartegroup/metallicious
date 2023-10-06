import warnings
import os
import MDAnalysis
import argparse
import parmed as pmd
import shutil
from tempfile import mkdtemp
from metallicious.log import logger
from metallicious.load_fingerprint import find_mapping_of_fingerprint_on_metal_and_its_surroundings
from metallicious.copy_topology_params import adjust_bonds, adjust_dihedrals, adjust_angles, adjust_impropers,\
    adjust_charge, adjust_pair_exclusions, update_pairs

warnings.filterwarnings('ignore')


class patcher():
    path = None
    tmpdir_path = None
    cage = None
    ligand = None
    metal_name = None
    metal_charge =None
    fp_topol = None
    fp_syst = None
    metal_indices =None
    topol_new = None # this is our new topology
    n_metals = None
    n_ligands = None
    name_of_binding_side=None
    fingerprint_style='full'

    output_topol = "cage.top"
    output_coords = "cage.gro"

    def __init__(self, cage_file=None, ligand_file=None, metal_name=None, metal_charge=None, name_of_binding_side=None,
                 cage_cgbind=None, smiles=None, arch_name=None, output_topol=None, output_coords=None):

        logger.info(f"Pathway of the library {os.path.dirname(__file__):s}")

        if output_topol is not None:
            self.output_topol = output_topol

        if output_coords is not None:
            self.output_coords = output_coords

        self.tmp_directory()
        logger.info('Current directory:' )
        logger.info(f'Created temporary directory: {self.tmpdir_path:s}')

    ''' TODO remove (2023/07/04)
    def from_cgbind(self): #TODO
        os.chdir(self.tmpdir_path)

    def from_smiles(self, smiles, arch_name, metal, metal_charge):
        os.chdir(self.tmpdir_path)

        self.create_cage_and_linker_cgbind(smiles, arch_name, metal, metal_charge)
        self.construct_cage(cage_file='cage.xyz', ligand_file='linker.top', metal_name=metal, metal_charge=metal_charge)
        #self.clean_up()
        os.chdir(self.path)

    def from_coords(self, cage_file, ligand_file, metal_name, metal_charge):
        shutil.copy(cage_file, f'{self.tmpdir_path:s}/{cage_file:s}')
        if ligand_file is not None:
            shutil.copy(ligand_file, f'{self.tmpdir_path:s}/{ligand_file:s}')
            
        os.chdir(self.tmpdir_path)
        self.construct_cage(cage_file=cage_file, ligand_file=ligand_file, metal_name=metal_name, metal_charge=metal_charge)
        #self.clean_up()
        os.chdir(self.path)
    '''

    def tmp_directory(self):
        self.path = os.getcwd()
        self.tmpdir_path = mkdtemp()

    def save(self,output_coords, output_topol, tmpdir_path):

        if output_coords.endswith('.gro'):
            shutil.copy(f'{self.cage_coord:s}', f'{output_coords:s}')
        else:
            coord = pmd.load_file(f'{tmpdir_path:s}/{self.cage_coord:s}')
            coord.save(f'{output_coords:s}', overwrite=True)

        self.topol_new.write(f"temp_topol.top")

        if output_topol.endswith('.top'):
            shutil.copy(f'temp_topol.top', f'{output_topol:s}')
        else:
            topol = pmd.load_file(f'temp_topol.top')
            topol.save(f'{output_topol:s}', overwrite=True)

        os.remove(f"temp_topol.top")

    def close(self):
        shutil.rmtree(self.tmpdir_path)
    ''' # TODO remove (2023/07/04)
    def construct_cage(self, cage_file=None, ligand_file=None, metal_name=None, metal_charge=None): #TODO remove
        
        The main function, which copies all the bonded paramters to the cage
        Sepearated for stages;
        1) Loads the cage
        2) Loads the fingerprint (it tries to make guess if now sure)
        3) Copies all the paramters

        :param cage_file:
        :param ligand_file:
        :param metal_name:
        :param metal_charge:
        :param name_of_binding_side:
        :return:
        

        self.metal_name = metal_name.title()
        self.metal_charge = metal_charge


        #self.load_cage(cage_file, ligand_file)
        self.prepare_new_topology(cage_file, self.metal_name, metal_charge=metal_charge, ligand_file=ligand_file)

        self.load_fingerprint_old(cage_file)

        for metal_index in self.metal_indices:
            mapping_fp_to_new, _ = find_mapping_of_fingerprint_on_metal_and_its_surroundings("cage.gro", metal_index, self.metal_name, self.fp_syst, cutoff=self.ligand_cutoff)

            self.adjust_charge(mapping_fp_to_new)

            self.adjust_bonds(mapping_fp_to_new)

            self.adjust_angles(mapping_fp_to_new)

            self.adjust_dihedrals(mapping_fp_to_new)

            self.adjust_impropers(mapping_fp_to_new)

        logger.info(f'Saving as {self.output_topol:s}')

        logger.info('Finished')
    '''


    def copy_site_topology_to_supramolecular(self, sites, cage_coord=None, cage_topol=None):
        '''
        The main function, which copies all the bonded paramters to the cage
        Sepearated for stages;
        1) Loads the cage
        2) Loads the fingerprint (it tries to make guess if now sure)
        3) Copies all the paramters

        :param cage_file:
        :param ligand_file:
        :param metal_name:
        :param metal_charge:
        :param name_of_binding_side:
        :return:
        '''

        #self.prepare_new_topology(cage_file, self.metal_name, metal_charge=metal_charge, ligand_file=ligand_file)

        self.cage_coord = cage_coord

        self.prepare_new_topology(cage_coord, cage_topol)

        for site in sites:

            mapping_fp_to_new, _ = find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_coord, site.index, site.metal_name, site.fp_syst, cutoff=site.ligand_cutoff, covalent_cutoff=site.covalent_cutoff)

            self.topol_new = adjust_charge(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_bonds(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_angles(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_dihedrals(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_impropers(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_pair_exclusions(self.topol_new, site.fp_topol, mapping_fp_to_new)


        # we need to remove pair exclusions, for atoms which are connected through angle containing metal
        self.topol_new = update_pairs(self.topol_new)

        logger.info('Finished')

    #def prepare_new_topology(self, cage_coord, cage_topol, metal_name, metal_charge, ligand_file=None):
    def prepare_new_topology(self, cage_coord, cage_topol):

        self.cage = MDAnalysis.Universe(cage_coord)

        # parametrize=False is needed to modify the parameters and save. However it does not load parameters from itp files
        # therfore we need to read it with parametrize=True and copy types of parameters
        self.topol_new = pmd.load_file(cage_topol, parametrize=False)

        topol_new2 = pmd.load_file(cage_topol)
        for idx, _ in enumerate(self.topol_new.atoms):
            self.topol_new.atoms[idx].type = topol_new2.atoms[idx].type
            self.topol_new.atoms[idx].epsilon = topol_new2.atoms[idx].epsilon
            self.topol_new.atoms[idx].sigma = topol_new2.atoms[idx].sigma
            self.topol_new.atoms[idx].rmin = topol_new2.atoms[idx].rmin
            self.topol_new.atoms[idx].charge = topol_new2.atoms[idx].charge

        for idx, _ in enumerate(self.topol_new.bonds):
            self.topol_new.bonds[idx].type = topol_new2.bonds[idx].type
        self.topol_new.bond_types = topol_new2.bond_types

        for idx, _ in enumerate(self.topol_new.angles):
            self.topol_new.angles[idx].type = topol_new2.angles[idx].type
        self.topol_new.angle_types = topol_new2.angle_types

        for idx, _ in enumerate(self.topol_new.dihedrals):
            self.topol_new.dihedrals[idx].type = topol_new2.dihedrals[idx].type
        self.topol_new.dihedral_types = topol_new2.dihedral_types

        for idx, _ in enumerate(self.topol_new.impropers):
            self.topol_new.impropers[idx].type = topol_new2.impropers[idx].type
        self.topol_new.improper_types = topol_new2.improper_types


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Structure pf the cage ")
    parser.add_argument("-linker_topol", help="Topology of the linker ") # TODO
    parser.add_argument("-linker_coords", help="Structure of the linker ")  # TODO
    parser.add_argument("-fingerprint", default=None, help="Structure of the linker ")  # TODO
    parser.add_argument("-fingerprint_style", default=None, help="Structure of the linker ")  # TODO

    parser.add_argument("-name", default='UNK', help="trajectory (traj_comp.xtc)")
    parser.add_argument("-smiles", default=None, help="Smiles of the linker")
    parser.add_argument("-arch_name", default=None, help="Architecture")
    parser.add_argument("-metal", default='Pd', help="Metal name")
    parser.add_argument("-metal_charge", default='2', help="Metal charge")

    parser.add_argument("-o", default='cage.gro', help="Output coordination file")
    parser.add_argument("-ot", default='cage.top', help="Output topology file")

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    cgbind2gmx = patcher()

    if args.fingerprint is not None:
        cgbind2gmx.name_of_binding_side = args.fingerprint


    if args.fingerprint_style is not None:
        cgbind2gmx.fingerprint_style = args.fingerprint_style

    if args.smiles is not None and args.arch_name is not None and args.metal is not None and args.metal_charge is not None:
        cgbind2gmx.from_smiles(args.smiles, args.arch_name, args.metal, int(args.metal_charge))
        cgbind2gmx.save(output_coords=args.o, output_topol=args.ot)
        cgbind2gmx.close()
    elif args.f is not None and args.linker_topol is not None and args.metal is not None and args.metal_charge is not None:
        cgbind2gmx.from_coords(args.f, args.linker_topol, args.metal, int(args.metal_charge))
        cgbind2gmx.save(output_coords=args.o, output_topol=args.ot)
        cgbind2gmx.close()
    else:
        print("One of the three values needs to be specified:")
        print('a) smiles, arch_name, metal, charge')
        print("b) f, linker_topol, metal, charge")


    #cgbind2gmx(cage_file="/media/piskorz/Oxford/Zincage/big_cage/cage.pdb", ligand_file='/media/piskorz/Oxford/Zincage/parameters/large_ligand/ZnL.top',
    #           metal_name = "ZN", metal_charge = 2, name_of_binding_side = "ZnN2")
