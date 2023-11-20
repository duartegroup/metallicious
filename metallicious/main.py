import warnings
import os
import MDAnalysis
import parmed as pmd
import shutil
from tempfile import mkdtemp
from metallicious.log import logger
from metallicious.load_fingerprint import find_mapping_of_fingerprint_on_metal_and_its_surroundings
from metallicious.copy_topology_params import adjust_bonds, adjust_dihedrals, adjust_angles, adjust_impropers,\
    adjust_charge, adjust_pair_exclusions, update_pairs, add_1_4_metal_pairs

warnings.filterwarnings('ignore')


class patcher():
    '''


    '''
    path = None
    tmpdir_path = None
    cage = None
    topol_new = None # this is our new topology

    output_topol = "cage.top"
    output_coords = "cage.gro"

    def __init__(self, output_topol=None, output_coords=None):
        None
        ''' This is redundant TODO (1/11/23)
        logger.info(f"Pathway of the library {os.path.dirname(__file__):s}")

        if output_topol is not None:
            self.output_topol = output_topol

        if output_coords is not None:
            self.output_coords = output_coords

        self.tmp_directory()
        
        logger.info('Current directory:' )
        logger.info(f'Created temporary directory: {self.tmpdir_path:s}')
        
    def tmp_directory(self):
        self.path = os.getcwd()
        self.tmpdir_path = mkdtemp()
        
    def close(self):
        shutil.rmtree(self.tmpdir_path)
        '''
    def save(self, output_coords, output_topol, tmpdir_path):
        if output_coords.endswith('.gro'):
            shutil.copy(f'{self.cage_coord:s}', f'{output_coords:s}')
        else:
            syst = MDAnalysis.Universe(f'{tmpdir_path:s}/{self.cage_coord:s}')
            syst.atoms.write(f'{output_coords:s}')

        self.topol_new.write(f"temp_topol.top")

        if output_topol.endswith('.top'):
            shutil.copy(f'temp_topol.top', f'{output_topol:s}')
        else:
            topol = pmd.load_file(f'temp_topol.top')
            topol.save(f'{output_topol:s}', overwrite=True)

        os.remove(f"temp_topol.top")



    def copy_site_topology_to_supramolecular(self, sites, cage_coord=None, cage_topol=None):
        '''
        The main function, which copies all the bonded paramters to the cage
        Sepearated for stages;
        1) Loads the cage
        2) Loads the fingerprint (it tries to make guess if now sure)
        3) Copies all the paramters

        :param sites:
        :param cage_coord:
        :param cage_topol:
        :return:
        '''

        self.cage_coord = cage_coord

        self.prepare_new_topology(cage_coord, cage_topol)


        for site in sites:
            mapping_fp_to_new, _ = find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_coord, site.index, site.metal_name, site.fp_syst, cutoff=site.ligand_cutoff, covalent_cutoff=site.covalent_cutoff, donors=site.donors)

            self.topol_new = adjust_charge(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_bonds(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_angles(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_dihedrals(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_impropers(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_pair_exclusions(self.topol_new, site.fp_topol, mapping_fp_to_new)

            metal_index = mapping_fp_to_new[0]
            self.topol_new = add_1_4_metal_pairs(self.topol_new, metal_index)


        # we need to remove pair exclusions, for atoms which are connected through angle containing metal
        self.topol_new = update_pairs(self.topol_new)
        logger.info('Finished')

    def prepare_new_topology(self, cage_coord, cage_topol):
        '''
        Copies existing topology to the topology which can be read and modified

        :param cage_coord: (string) filename of the coordination file of the structure
        :param cage_topol: (string) filename of the topology (readable by parmed) file of the structure
        :return:
        '''

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

