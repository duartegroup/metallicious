from metallicious.extract_metal_site import extract_metal_structure, find_metal_indices
from metallicious.seminario import single_seminario, check_if_orca_available
from metallicious.charges import calculate_charges2
from metallicious.copy_topology_params import copy_bonds, copy_angles, copy_dihedrals, copy_impropers, \
    copy_pair_exclusions, update_pairs
from metallicious.load_fingerprint import guess_fingerprint, load_fp_from_file
from metallicious.data import vdw_data
from metallicious.prepare_initial_topology import prepare_initial_topology
from metallicious.log import logger
from metallicious.patcher import patcher
from metallicious.utils import new_directory
from metallicious.asserts import compare_topology_and_coords, check_if_parametrization_modules_available

import numpy as np
import MDAnalysis
from MDAnalysis.lib.distances import distance_array

import shutil
import os

class supramolecular_structure:
    '''
    The main structure holidng all the information about the intput and communicating between different classes

    '''
    def __init__(self, filename, metal_charge_mult=None, metal_charges=None, LJ_type=None, topol=None,
                 keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], improper_metal=False,
                 donors=['N', 'S', 'O'],
                 library_path=f'{os.path.dirname(__file__):s}/library/', ff='gaff', search_library=True,
                 fingerprint_guess_list=None, truncation_scheme=None, covalent_cutoff=3, rmsd_cutoff=2):
        '''
        Initialize the class

        :param filename: (str) name of the coordination file
        :param metal_charge_mult:  (dict) the names charges, and multiplicity of the metals in format
                                        {metal_name: (metal_charges, multiplicity)}
        :param metal_charges: (dict) the names and charges of metals in the input structure in format:
                                          {metal_name1: metal_charges1, metal1_name2: metal_charge2}
        :param LJ_type: (str) name of LJ dataset used for metal paramters
        :param topol: (str) path to topology (optional)
        :param keywords: list(str) the keywords for QM calculations
        :param improper_metal: (bool) if True it will parametrize the improper dihedral involving metal
        :param donors: (list(str)) list of atom elements with which metal forms bond
        :param library_path: (str) directory of template library, be default where the script is
        :param ff: (str) parametrization protocol for small organic molecules (only gaff available)
        :param search_library: (bool) if True, metallicious searches templates in template library,
                    if False, it will parametrize template
        :param fingerprint_guess_list: (list(str)) list of templates to check
        :param truncation_scheme: (str) name of the truncation scheme
        :param covalent_cutoff: (float) if metal-atoms smaller then cutoff it creates bonds ligand with metal
        :param rmsd_cutoff: (float) cutoff for the RMSD acceptance of the template
        '''

        logger.info(f"Library with templates is located: {library_path:}")

        self.unique_sites = []
        self.sites = []
        self.filename = filename
        self.library = library_path

        self.autode_keywords = keywords
        self.donors = donors
        self.improper_metal = improper_metal
        self.ff = ff

        self.covalent_cutoff = covalent_cutoff
        self.closest_neighbhors = 3

        self.rmsd_cutoff = rmsd_cutoff

        #if topol is not None:
        #self.topol = self.make_metals_first(topol)
        self.topol = topol

        self.allow_new_templates = True

        self.metal_charge_dict = {}
        self.metal_mult_dict = None

        if metal_charge_mult is not None:
            self.metal_charge_dict = {}
            self.metal_mult_dict = {}
            for key in metal_charge_mult:
                self.metal_charge_dict[key] = metal_charge_mult[key][0]
                self.metal_mult_dict[key] = metal_charge_mult[key][1]

            self.metal_names = list(metal_charge_mult.keys())
        elif metal_charges is not None:
            self.metal_charge_dict = metal_charges
            self.metal_names = list(metal_charges.keys())
            self.allow_new_templates = False
        else:
            raise ValueError("Not correct format of metal_charge_mult/metal_charges")

        if LJ_type == 'custom' and topol is not None:
            self.vdw_type = 'custom'

        elif LJ_type is None:
            for LJ_type in vdw_data:
                present = [(metal_name in vdw_data[LJ_type] or f"{metal_name:}{self.metal_charge_dict[metal_name]:}" in vdw_data[LJ_type]) for metal_name in self.metal_names]
                if sum(present) == len(present):
                    self.vdw_type = LJ_type
                    logger.info(
                        f"vdw_type not selected, will use first available for selected metals: {LJ_type:}")
        else:
            for metal_name in self.metal_names:
                if metal_name not in vdw_data[LJ_type] and f"{metal_name:}{self.metal_charge_dict[metal_name]:}" not in \
                        vdw_data[LJ_type]:

                    metal_available_in = []
                    for LJ_type in vdw_data:
                        if (metal_name in vdw_data[LJ_type] or f"{metal_name:}{self.metal_charge_dict[metal_name]:}" in vdw_data[LJ_type]):
                            metal_available_in.append(LJ_type)

                    raise ValueError(f"Metal ({metal_name:}) unavailable in selected LJ library, but it seems it is present in {metal_available_in:}")
            self.vdw_type = LJ_type

        self.fingerprint_guess_list = fingerprint_guess_list
        self.truncation_scheme = truncation_scheme
        self.search_library = search_library
        self.library_path = f"{os.path.dirname(__file__):s}/library"

        self.find_metal_sites()

        self.path = os.getcwd()
        # This is placeholder for temporary direction, but ORCA often breaks, and it is easier just to restart a job
        self.tmpdir_path = '.'  # mkdtemp()
        os.chdir(self.tmpdir_path)


    def find_metal_sites(self):
        syst = MDAnalysis.Universe(self.filename)
        for name in self.metal_names:
            indices, _ = find_metal_indices(syst, name)
            for index in indices:
                site = metal_site(name.title(), self.metal_charge_dict[name], index, fp_style=self.truncation_scheme, covalent_cutoff=self.covalent_cutoff, donors=self.donors)
                self.sites.append(site)
        self.assign_fingerprints()


    def assign_fingerprints(self):
        additional_fp_coords = {}
        for unique_site in self.unique_sites:
            if unique_site.fp_coord_file is not None:
                suffix = '_' + unique_site.directory.split('_')[1]
                additional_fp_coords[unique_site.name+suffix] = unique_site.fp_coord_file

        logger.info(f"Fingerprint to choose from: {additional_fp_coords:}")
        for site in self.sites:
            guessed = guess_fingerprint(self.filename, site.index, metal_name=site.metal_name,
                                        metal_charge=site.metal_charge,
                                        fingerprint_guess_list=self.fingerprint_guess_list,
                                        m_m_cutoff=10, vdw_type=self.vdw_type, library_path=self.library_path,
                                        search_library=self.search_library,
                                        additional_fp_files=additional_fp_coords, fp_style=self.truncation_scheme,
                                        donors=self.donors, rmsd_cutoff=self.rmsd_cutoff)

            if guessed is not False:  # do not change to True...
                site.fp_coord_file = f"{guessed:}.pdb"
                site.fp_topol_file = f"{guessed:}.top"
                site.load_fingerprint()
                site.set_cutoff()
            else:
                logger.info("Template for this site not found")

    def extract_unique_metal_sites(self):
        logger.info(f"Extracting")
        unique_sites = []

        # extract metal sites:
        for metal_name in self.metal_names:
            site_lists = extract_metal_structure(self.filename, self.topol, metal_name, output=f"site_{metal_name:}",
                                                 all_metal_names=self.metal_names, covalent_cutoff=self.covalent_cutoff, donors=self.donors, closest_neighbhors = self.closest_neighbhors)

            for site_list in site_lists:
                site_list[1] = self.metal_charge_dict[metal_name]  # we change the charge
                if self.metal_mult_dict is not None:
                    site_list[2] = self.metal_mult_dict[metal_name]  # we change the multiplicity

                site_list += [self.autode_keywords, self.improper_metal, self.donors, self.vdw_type, self.rmsd_cutoff]
                unique_sites += [new_metal_site(*site_list)]

        self.unique_sites = unique_sites

    def check_if_parameters_available(self):
        #if len([1 for site in self.sites if site.fp_topol_file is not None]) == len(self.sites):
        if len([1 for site in self.sites if site.fp_topol_file is not None]) == len(self.sites):
            return True
        else:
            return False

    def prepare_initial_topology(self, coord_filename='noncovalent_complex.pdb',
                                 topol_filename='noncovalent_complex.top', method='gaff', homoleptic_ligand_topol=None,
                                 subdir='init_topol'):
        here = os.getcwd()
        new_directory(subdir)

        old_filename = self.filename
        new_filename = self.filename[self.filename.rfind('/')+1:]

        try:
            shutil.copyfile(old_filename, f'{subdir}/{new_filename}')
        except shutil.SameFileError:
            pass

        if homoleptic_ligand_topol is not None:
            shutil.copyfile(homoleptic_ligand_topol, f'{subdir}/{homoleptic_ligand_topol}')
        os.chdir(subdir)

        if method == 'gaff' or homoleptic_ligand_topol is not None:
            metal_indicies = prepare_initial_topology(new_filename, self.metal_names, self.sites[0].metal_charge,
                                                      coord_filename, topol_filename, self.vdw_type,
                                                      ligand_topol=homoleptic_ligand_topol)
            self.filename = f'{subdir}/{coord_filename}'
            self.topol = f'{subdir}/{topol_filename}'

            for metal_index, site in zip(metal_indicies, self.sites):
                site.index = metal_index
        else:
            raise ValueError("Only GAFF supported")
        os.chdir(here)

    def parametrize(self, out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=False):
        if prepare_initial_topology==True:
            self.prepare_initial_topology()

        if self.topol is None:
            raise Exception("Topology file not specified, please provide topology, or use prepare_initial_topology=True")

        if self.check_if_parameters_available() is False:
            logger.info("[ ] Extracting the structure")
            self.parametrize_metal_sites()
            self.assign_fingerprints()

        logger.info(self.summary())
        logger.info("The templates are available!")

        # Check if number of atoms agree in both files
        compare_topology_and_coords(self.topol, self.filename)

        parameter_copier = patcher()
        parameter_copier.copy_site_topology_to_supramolecular(self.sites, cage_coord=self.filename,
                                                              cage_topol=self.topol)

        parameter_copier.save(f'{self.path:s}/{out_coord:s}', f'{self.path:s}/{out_topol:s}',
                              tmpdir_path=self.tmpdir_path)
        logger.info("[+] Finished!")
        return True

    def parametrize_metal_sites(self):
        if len(self.unique_sites) == 0:
            self.extract_unique_metal_sites()

        if len(self.unique_sites)>0 and self.allow_new_templates is True:
            check_if_orca_available()

        for site in self.unique_sites:
            #if site.check_library() is False:
            if site.fp_topol_file is None:
                if self.allow_new_templates is True:
                    check_if_parametrization_modules_available()

                    site.parametrize()
                    self.add_site_to_library(site)
                else:
                    raise Exception(
                        "Template not found (try to (a) parametrize it (specify multiplicity) or (b) truncate template)")


    def add_site_to_library(self, site):
        # adding to the library
        if os.path.isfile(site.fp_topol_file) and os.path.isfile(site.fp_coord_file):
            file_idx = 0
            while True:
                if os.path.isfile(f"{self.library_path:}/{site.name}_{file_idx}.top"):
                    file_idx += 1
                else:
                    break

            if self.vdw_type != 'custom':  # if it custom we don't want it
                logger.info(f"[+] Saving as {self.library_path:}/{site.name}_{file_idx}.top")
                shutil.copyfile(site.fp_topol_file, f"{self.library_path:}/{site.name}_{file_idx}.top")
                shutil.copyfile(site.fp_coord_file, f"{self.library_path:}/{site.name}_{file_idx}.pdb")
                file_idx += 1

    def summary(self):
        string = "Sites:\n"
        for site in self.sites:
            string += site._print() + '\n'

        string +=  "\nUnique sites:"
        for site in self.unique_sites:
            string += f"{site.metal_name}({site.metal_charge}+)"

class metal_site():
    '''
    This class stores information about the metal binding site part of the input structure
    '''

    def __init__(self, metal_name, metal_charge, index, fp_topol=None, fp_coord=None, fp_style=None, covalent_cutoff=3.0, donors = None):
        '''
        Information needed later for patcher to modify the input force-field parameters

        :param metal_name: (names) metal names
        :param metal_charge: (int) formal of charge
        :param index: (int) index of metal
        :param fp_topol: (parmed.topology) topology of template
        :param fp_coord: (MDAnalysis.Universe) coordinates of the template
        :param fp_style: (str) truncation scheme for template
        :param covalent_cutoff: (float) cut-off for metal-ligand interactions
        :param donors: (list(str)) list of elements used as donor-acceptor
        '''

        self.metal_name = metal_name
        self.metal_charge = metal_charge
        self.index = index
        self.fp_topol_file = fp_topol
        self.fp_coord_file = fp_coord
        self.fp_style = fp_style
        self.covalent_cutoff= covalent_cutoff
        self.ignore_truncation_warning = True
        self.donors = donors

    def _print(self):
        if self.fp_coord_file is not None:
            return f"<{self.index}: {self.metal_name}({self.metal_charge}+) {self.fp_coord_file.split('/')[-1]}>"
        else:
            return f"<{self.index}: {self.metal_name}({self.metal_charge}+) None>"

    def __str__(self):
        return self._print()

    def __repr__(self):
        return self._print()

    def load_fingerprint(self):
        '''
        Loads files from the library into the class, if the name is not known, it will try to guess the fingerprint

        :param name_of_binding_side:
        :return:
        '''

        self.fp_topol, self.fp_syst = load_fp_from_file(self.fp_coord_file, self.fp_topol_file, self.fp_style, self.ignore_truncation_warning)
        return True

    def set_cutoff(self):
        metal_position = self.fp_syst.atoms[0].position
        nometal_position = self.fp_syst.atoms[1:].positions
        self.ligand_cutoff = np.max(distance_array(metal_position, nometal_position)) + 2.0

        return True


class new_metal_site():
    def __init__(self, metal_name, metal_charge, mult, directory, ligand_names=None, unique_ligands_pattern=None,
                 link_atoms=None, additional_atoms=None, starting_index=None, indecies=None, ligand_charges=None,
                 ligand_smiles=None, topol=None,
                 keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], improper_metal=True,
                 donors=['N', 'S', 'O'], vdw_type='uff', rmsd_cutoff=2):
        '''
        The class stores and processes information needed for new template parametrization

        :param metal_name: (str) metal name
        :param metal_charge: (int) formal charge
        :param mult: (int) multiplicity
        :param directory: (str) directory in which parameterization will be performed
        :param ligand_names: (list(str)) list of ligand names
        :param unique_ligands_pattern: (list(int)) list indicated which ligands are the same (e.g. [0,1,1])
        :param link_atoms: (list(int)) indices of link atoms (unsaturated atoms)
        :param additional_atoms: (list(int)) indices of additional hydrogen's in saturate template, which will be removed
        :param starting_index: (list(int)) starting indices for new ligands
        :param indecies: list(int) indices of template's atoms in the input structure
        :param ligand_charges: (list(int)) list of ligand charges
        :param ligand_smiles: (list(str)) list of ligand smiles
        :param topol: (str) (Non-bonded) force-field parameters of the metal site (which will be changed)
        :param keywords: list(str) keywords for QM cancellations
        :param improper_metal: (bool) if True improper dihedral involving metal will be parametrized
        :param donors: (list(str)) list of element names with which metal are connected
        :param vdw_type: (str) name of the dataset used for metal Lennard-Jones parameters
        :param rmsd_cutoff: (float) cutoff for the RMSD acceptance of the template
        '''

        self.unique_ligands_pattern = unique_ligands_pattern
        self.link_atoms = link_atoms
        self.additional_atoms = additional_atoms
        self.starting_index = starting_index
        self.indecies = indecies
        self.metal_name = metal_name
        self.metal_radius = None

        self.metal_charge = metal_charge
        self.mult = mult
        self.ligand_charges = ligand_charges
        # smiles are not used they are just for user if manual correction needs to be done
        # it is easier to look at the specific ligand to change its charge etc.
        self.ligand_smiles = ligand_smiles

        self.bonds = None
        self.angles = None
        self.dihedrals = None

        self.directory = directory
        self.filename = 'saturated_template.xyz'
        self.ligand_names = ligand_names
        self.topol = None

        self.autode_keywords = keywords
        self.donors = donors
        self.improper_metal = improper_metal
        self.vdw_type = vdw_type

        self.fp_name = None
        self.fp_topol_file = None
        self.fp_coord_file = None

        self.name = f"{metal_name:}_{metal_charge:}_{vdw_type:}"
        self.topol = topol

        self.vibrational_scaling = None

        self.rmsd_cutoff = rmsd_cutoff

        if vdw_type is not None:
            if metal_name in vdw_data[vdw_type]:
                vdw_entry = metal_name
            elif f"{metal_name:}{metal_charge:}" in vdw_data[vdw_type]:
                vdw_entry = f"{metal_name:}{metal_charge:}"
            eps, r2min = vdw_data[vdw_type][vdw_entry]
            self.metal_radius = r2min # *2 radius

            # Change metal type:
            self.topol[0].atom_type.rmin = r2min
            self.topol[0].atom_type.epsilon = eps
        else:
            if self.topol is not None:
                self.metal_radius = self.read_radius_from_topol()

    def read_radius_from_topol(self):
        self.metal_radius = self.topol[0].rmin

    def _print(self):
        return f"<{self.metal_name}({self.metal_charge}+) {self.name}>"

    def __str__(self):
        return self._print()

    def __repr__(self):
        return self._print()

    def seminario(self):
        '''
        Find the distances and force constants of metal binding site using Cole's Seminario method
        :return:
        '''
        site_charge = self.metal_charge + np.sum(self.ligand_charges)

        self.bonds, self.angles, self.dihedrals, self.impropers, self.pairs, self.filename = single_seminario(
            self.filename, site_charge,
            self.metal_name,
            self.starting_index, self.indecies,
            self.unique_ligands_pattern,
            keywords=self.autode_keywords,
            mult=self.mult,
            improper_metal=self.improper_metal,
            donors=self.donors,
            atoms_to_remove=self.additional_atoms,
            vibrational_scaling=self.vibrational_scaling)

        self.topol = copy_bonds(self.topol, self.bonds, self.metal_name)
        self.topol = copy_angles(self.topol, self.angles, self.metal_name)

        if len(self.dihedrals) > 0:
            self.topol = copy_dihedrals(self.topol, self.dihedrals, self.metal_name)
        if len(self.impropers) > 0:
            self.topol = copy_impropers(self.topol, self.impropers, self.metal_name)

        if len(self.pairs) > 0:
            self.topol = copy_pair_exclusions(self.topol, self.pairs)



        self.topol = update_pairs(self.topol)


    def partial_charge(self):
        '''
        Calculate the charges of the metal site using restrained electrostatic poential (RESP)
        :return:
        '''

        partial_charges = calculate_charges2(self.metal_name, self.metal_charge, self.directory,
                                             self.unique_ligands_pattern, self.ligand_charges, self.link_atoms,
                                             self.additional_atoms,
                                             self.starting_index,
                                             metal_radius=self.metal_radius, mult=self.mult)

        # Removing additional atoms:
        partial_charges = [partial_charge for idx, partial_charge in enumerate(partial_charges) if
                           idx not in self.additional_atoms]

        logger.info("\t[ ] Copying charges")
        for idx, atom in enumerate(self.topol.atoms):
            atom.charge = partial_charges[idx]
        logger.info("\t[+] Charges calculated !")


    def reduce_to_template(self, out_topol='template.top', out_coord='template.pdb'):
        '''
        Remove the extra atoms(hydrogens) from saturated template to create template

        :param out_topol: (str) output force-field parameters file
        :param out_coord: (str) ouput coordination file
        :return:
        '''
        # self.topol.write(f"old_new_topol.top")

        # self.topol.strip(f"@{','.join(list(map(str, np.array(self.additional_atoms) + 1))):s}")
        self.topol.save(out_topol, overwrite=True)

        new_cage = MDAnalysis.Universe(self.filename)
        new_cage.atoms.write(f"old_new_template.pdb")

        logger.info(f"removing additional atoms {self.additional_atoms:}")

        if len(self.additional_atoms) > 0:
            template = new_cage.select_atoms(f"not index {' '.join(list(map(str, np.array(self.additional_atoms)))):s}")
        else:
            template = new_cage.atoms

        template.write(out_coord)

        self.fp_topol_file = f"{self.directory}/{out_topol}"
        self.fp_coord_file = f"{self.directory}/{out_coord}"


        logger.info(f"Topology file: {self.fp_topol_file:}")
        logger.info(f"Coordination file: {self.fp_coord_file:}")

    def check_library(self):
        '''
        Guess the template of the metal site
        :return:
        '''
        return guess_fingerprint(self.directory + "/saturated_template.xyz", 0, metal_name=self.metal_name, metal_charge=self.metal_charge, vdw_type=self.vdw_type, rmsd_cutoff=self.rmsd_cutoff)

    def parametrize(self):
        '''
        Parametrizes the new template: runs seminarion, partial charge calculation and reduction to the template
        :return:
        '''
        here = os.getcwd()
        os.chdir(self.directory)

        self.seminario()
        self.partial_charge()
        self.reduce_to_template()
        os.chdir(here)
        logger.info("Site parametrized successfully!")
