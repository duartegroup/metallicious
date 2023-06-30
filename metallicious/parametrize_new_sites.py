import argparse
import shutil

# TODO this should be class not INFO file
from extract_metal_site import extract_metal_structure, read_info_file, find_metal_indices
from seminario import single_seminario
from charges import calculate_charges2
from initial_site import create_initial_topol2
from copy_topology_params import copy_bonds, copy_angles, copy_dihedrals
from load_fingerprint import guess_fingerprint, load_fp_from_file
from data import vdw_data
from prepare_initial_topology import prepare_initial_topology

from tempfile import mkdtemp



import parmed as pmd
import numpy as np
import MDAnalysis
from MDAnalysis.lib.distances import distance_array

import os
from main import cgbind2pmd

from log import logger


class supramolecular_structure:
    def __init__(self, filename, metal_charge_mult=None, metal_charges=None, vdw_type=None, topol=None,
                 keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], improper_metal=True,
                 donors=['N', 'S', 'O'],
                 library_path=f'{os.path.dirname(__file__):s}/library/', ff='gaff', search_library=True,
                 fingerprint_guess_list=None, fingerprint_style = 'full', new_template=True):

        self.unique_sites = []
        self.sites = []
        self.filename = filename
        self.library = library_path

        self.autode_keywords = keywords
        self.donors = donors
        self.improper_metal = improper_metal
        self.ff = ff

        self.topol = topol

        if metal_charge_mult is None:
            print("Guessing metals not implemented ")  # TODO
        else:
            self.metal_charge_dict = metal_charge_mult
            self.metal_names = list(metal_charge_mult.keys())

        if vdw_type is None and topol is not None:
            self.vdw_type = 'custom'
            None
            # use vdw from topol


        elif vdw_type is None:
            for vdw_type in vdw_data:
                present = [metal_name in vdw_data[vdw_type] for metal_name in self.metal_names]
                if sum(present) == len(present):
                    self.vdw_type = vdw_type
                    print(f"vdw_type not selected, will use first available for selected metals: {vdw_type:}")
                    break
        #elif vdw_type == 'force':
        #    self.vdw_type = None
        else:
            for metal_name in self.metal_names:
                if metal_name not in vdw_data[vdw_type]:
                    raise

            self.vdw_type = vdw_type

        self.fingerprint_guess_list = fingerprint_guess_list
        self.fingerprint_style = 'full'
        self.search_library = search_library
        self.library_path = f"{os.path.dirname(__file__):s}/library"

        self.find_metal_sites()

        # copy main file and create temp directory

        self.path = os.getcwd()


        self.tmpdir_path = '..'  #mkdtemp()
        #shutil.copy(self.filename, f'{self.tmpdir_path:}')
        os.chdir(self.tmpdir_path)

        if self.topol is None:
            print(f"Please provide topology, or use .prepare_initial_topology()")

    def find_metal_sites(self):
        syst = MDAnalysis.Universe(self.filename)
        for name in self.metal_names:

            #TODO

            indices, _ = find_metal_indices(syst, name)

            #metal_atoms = syst.select_atoms(f'name {name.title():}* {name.upper():}* {name.lower():}*')
            for index in indices:
                site = metal_site(name.title(), self.metal_charge_dict[name][0], index)
                self.sites.append(site)

        self.assign_fingerprints()

    def assign_fingerprints(self):
        additional_fp_coords = {}
        for unique_site in self.unique_sites:
            if unique_site.fp_coord_file is not None:
                additional_fp_coords[unique_site.name] = unique_site.fp_coord_file

        for site in self.sites:
            guessed = guess_fingerprint(self.filename, site.index, metal_name=site.metal_name,
                                        metal_charge=site.metal_charge,
                                        fingerprint_guess_list=self.fingerprint_guess_list,
                                        m_m_cutoff=10, vdw_type=self.vdw_type, library_path=self.library_path,
                                        search_library=self.search_library,
                                        additional_fp_files=additional_fp_coords) #TODO what about fp style (?)

            if guessed is not False:
                site.fp_coord_file = f"{self.library_path:s}/{guessed:}.pdb" #TODO
                site.fp_topol_file = f"{self.library_path:s}/{guessed:}.top"
                site.load_fingerprint()
                site.set_cutoff()

    def extract_unique_metal_sites(self):
        logger.info(f"Extracting")

        # Resetting Info file and sites
        if os.path.isfile("INFO.dat"):  # TODO is info file still needed?
            os.remove('INFO.dat')

        unique_sites = []

        if self.topol is None:
            print("please provide topology or use .prepare_initial_topol()")
            return False

        # extract metal sites:
        for metal_name in self.metal_names:
            site_lists = extract_metal_structure(self.filename, self.topol, metal_name, output=f"site_{metal_name:}",
                                                 all_metal_names=self.metal_names)

            for site_list in site_lists:
                site_list[1] = self.metal_charge_dict[metal_name][0]  # we change the charge
                site_list[2] = self.metal_charge_dict[metal_name][1]  # we change the multiplicity
                site_list += [self.autode_keywords, self.improper_metal, self.donors, self.vdw_type]
                unique_sites += [new_metal_site(*site_list)]

        self.unique_sites = unique_sites

    def read_info_file(self, filename='INFO.dat'):
        self.unique_sites = read_info_file(filename)

    def check_if_parameters_available(self):
        if len([1 for site in self.sites if site.fp_topol_file is not None]) == len(self.sites):
            return True
        else:
            return False

    def prepare_initial_topology(self, coord_filename= 'noncovalent_complex.pdb', topol_filename = 'noncovalent_complex.top', method='gaff', homoleptic_ligand_topol=None):
        if method=='gaff':
            metal_indicies = prepare_initial_topology(self.filename, self.metal_names, self.sites[0].metal_charge,
                                                      coord_filename, topol_filename, ligand_topol=homoleptic_ligand_topol)
            self.filename = coord_filename
            self.topol = topol_filename

            for metal_index, site in zip(metal_indicies, self.sites):
                site.index = metal_index
        else:
            print("Only gaff suported")

    def parametrize(self, out_coord='out.pdb', out_topol='out.top', parametrize_metal_sites=True):
        print("Parametrizing")
        if self.check_if_parameters_available() is False:
            print("[ ] Extracting the structure")
            self.parametrize_metal_sites()
            self.assign_fingerprints()

        self.summary()
        print("Available!")


        if self.topol is None:
            print(f"Please provide topology, or use .prepare_initial_topology()")
            return False


        else:
            print("Checking topology") # TODO Metals needs to be at the beggining of the file
            '''
            topol = pmd.load_file(self.topol)
            #topol.write(f"complex.top", [list(range(len(topol.split())))])
            new_top = topol.split()
            (new_top[1][0]+new_top[0][0]).write("complex.top", [[1,0]])
            self.topol = "complex.top"
            '''
            # TODO do we need metal be first ?




        print(self.topol)


        if self.topol is None:
            print("No topology")
            raise

        # check if topol and coord have the same names of atoms TODO

        parameter_copier = cgbind2pmd()
        parameter_copier.copy_site_topology_to_supramolecular(self.sites, cage_coord=self.filename, cage_topol=self.topol)
        parameter_copier.save(f'{self.path:s}/{out_coord:s}', f'{self.path:s}/{out_topol:s}', tmpdir_path=self.tmpdir_path)
        print("Finished!")
        return True


    def parametrize_metal_sites(self):
        #TODO check if ORCA available !!
        if len(self.unique_sites) == 0:
            self.extract_unique_metal_sites()

        for site in self.unique_sites:
            if site.check_library() is False:
                site.parametrize()
                self.add_site_to_library(site)


    def add_site_to_library(self, site):
        # adding to the library
        if os.path.isfile(site.fp_topol_file) and os.path.isfile(site.fp_coord_file):
            file_idx = 0

            while True:
                if os.path.isfile(f"{self.library_path:}/{site.name}_{file_idx}.top"):
                    file_idx+=1
                else:
                    break

            if self.vdw_type != 'custom': # if it custom we don't want it
                print(f"Saving as {self.library_path:}/{site.name}_{file_idx}.top")
                shutil.copyfile(site.fp_topol_file, f"{self.library_path:}/{site.name}_{file_idx}.top")
                shutil.copyfile(site.fp_coord_file, f"{self.library_path:}/{site.name}_{file_idx}.pdb")
                file_idx += 1

    def summary(self):
        print("Sites:")
        for site in self.sites:
            print(site)

        print("\nUnique sites:")
        for site in self.unique_sites:
            print(f"{site.metal_name}({site.metal_charge}+)")



class metal_site():
    def __init__(self, metal_name, metal_charge, index, fp_topol=None, fp_coord=None, fingerprint_style='full'):
        self.metal_name = metal_name
        self.metal_charge = metal_charge
        self.index = index
        self.fp_topol_file = fp_topol
        self.fp_coord_file = fp_coord
        self.fp_style = fingerprint_style # TODO chang all names "finerprints" to "fp"

    def _print(self):
        if self.fp_coord_file is not None:
            return f"{self.index}: {self.metal_name}({self.metal_charge}+) {self.fp_coord_file.split('/')[-1]}"
        else:
            return f"{self.index}: {self.metal_name}({self.metal_charge}+) None"

    def __str__(self):
        return self._print()
    def __repr__(self):
        return self._print()

    def load_fingerprint(self, fp_style='full'):
        '''
        Loads files from the library into the class, if the name is not known, it will try to guess the fingerprint

        :param name_of_binding_side:
        :return:
        '''

        # Inputd
        self.fp_topol, self.fp_syst = load_fp_from_file(self.fp_coord_file, self.fp_topol_file, fp_style)


        return True

    def set_cutoff(self):
        metal_position = self.fp_syst.atoms[0].position
        nometal_position = self.fp_syst.atoms[1:].positions

        self.ligand_cutoff = np.max(distance_array(metal_position, nometal_position)) + 2.0

        return True


        '''
        # In past I used this one, which should be working, but I think, I could solved the issue with ligands on the way. So, checked it and delete after testes:
        
        if self.n_metals > 1:
            mm_distances = distance_array(self.cage.atoms[self.metal_indices].positions, self.cage.atoms[self.metal_indices].positions)
            self.m_m_cutoff = np.min(mm_distances[mm_distances>0.1])-0.1 #this acctually does not make it better TODO
        else:
            self.m_m_cutoff = 1e10
            #self.m_m_cutoff=9            
        '''

class new_metal_site():
    def __init__(self, metal_name, metal_charge, mult, directory, ligand_names=None, unique_ligands_pattern=None,
                 link_atoms=None, additional_atoms=None, starting_index=None, indecies=None, ligand_charges=None,
                 ligand_smiles=None, topol=None,
                 keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'], improper_metal=True,
                 donors=['N', 'S', 'O'], vdw_type=None):
        
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
        self.ligand_smiles = ligand_smiles  # smiles are not really used, they are just for user if hand correction needs to be done

        self.bonds = None
        self.angles = None
        self.dihedrals = None

        self.directory = directory
        self.filename = 'site.xyz' # TODO maybe this needs to be also passed
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

        if vdw_type is not None:
            self.metal_radius =  vdw_data[vdw_type][metal_name][1]
        else:
            if self.topol is not None:
                self.metal_radius = self.read_radius_from_topol()

    def read_radius_from_topol(self): #TODO
        self.metal_radius = self.topol[0].rmin

    def _print(self):
        return f"{self.metal_name}({self.metal_charge}+) {self.name}"

    def __str__(self):
        return self._print()
    def __repr__(self):
        return self._print()

    def create_initial_topol(self): # TODO depracticed (22/06/2023)
        self.topol = create_initial_topol2(self.metal_name, self.metal_charge, self.unique_ligands_pattern,
                                           self.vdw_type)

    def seminario(self):
        site_charge = self.metal_charge + np.sum(self.ligand_charges)

        self.bonds, self.angles, self.dihedrals, self.filename = single_seminario(self.filename, site_charge,
                                                                                   self.metal_name,
                                                                                   self.starting_index, self.indecies,
                                                                                   self.unique_ligands_pattern,
                                                                                   keywords=self.autode_keywords,
                                                                                   mult=self.mult,
                                                                                   improper_metal=self.improper_metal,
                                                                                   donors=self.donors,
                                                                                  atoms_to_remove = self.additional_atoms)

        self.topol = copy_bonds(self.topol, self.bonds, self.metal_name)
        self.topol = copy_angles(self.topol, self.angles, self.metal_name)
        self.topol = copy_dihedrals(self.topol, self.dihedrals, self.metal_name)

    def partial_charge(self):
        # calculate_charges2(metal_name, metal_charge, filename, unique_ligands_pattern, link_atoms, additional_atoms,
        #                        starting_index, vdw_data_name, mult=1):
        partial_charges = calculate_charges2(self.metal_name, self.metal_charge, self.directory,
                                             self.unique_ligands_pattern, self.ligand_charges, self.link_atoms,
                                             self.additional_atoms,
                                             self.starting_index,
                                             metal_radius=self.metal_radius, mult=self.mult) # TODO vdw should be not specified,

        #Removing additional atoms:
        partial_charges = [partial_charge for idx, partial_charge in enumerate(partial_charges) if idx not in self.additional_atoms]
        
        print("\t[ ] Copying charges")
        for idx, atom in enumerate(self.topol.atoms):
            atom.charge = partial_charges[idx]
        print("\t[+] Charges calculated !")

    def reduce_to_template(self, out_topol='template.top', out_coord='template.pdb'):
        #self.topol.write(f"old_new_topol.top")

        #self.topol.strip(f"@{','.join(list(map(str, np.array(self.additional_atoms) + 1))):s}")
        self.topol.save(out_topol, overwrite=True)

        new_cage = MDAnalysis.Universe(self.filename)
        #new_cage.atoms.write(f"old_new_template.pdb")

        new_cage.select_atoms(f"not index {' '.join(list(map(str, np.array(self.additional_atoms)))):s}").write(
            out_coord)  # saving should be speartate TODO

        self.fp_topol_file = f"{self.directory}/{out_topol}"
        self.fp_coord_file = f"{self.directory}/{out_coord}"

    def check_library(self):
        return guess_fingerprint(self.directory + "/site.xyz", 0, metal_name=self.metal_name)  # TODO library

    def parametrize(self):
        here = os.getcwd()
        os.chdir(self.directory)

        #self.create_initial_topol()
        self.seminario()
        self.partial_charge()
        self.reduce_to_template()
        os.chdir(here)


        print("Site parametrized successfully!")


def parametrize(filename, metal_name, metal_charge, keywords=['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'],
                vdw_type="uff", mult=1, improper_metal=True, donors=['N', 'S', 'O']):
    """
    Combines all other scripts to extract, and parameterize metal site(s)
    :param filename:
    :param metal_name:
    :param metal_charge:
    :return:
    """

    if os.path.isfile("INFO.dat"):
        os.remove('INFO.dat')

    print("[ ] Extracting the structure")
    extract_metal_structure(filename, metal_name, output="site")
    print("[ ] Create initial topology")
    topols = create_initial_topol(metal_name, metal_charge, vdw_data_name=vdw_type)
    n_sites = len(topols)

    topols = []
    for n_site in range(n_sites):
        topols.append(pmd.load_file(f"topol{n_site:d}.top"))
    print(topols)

    print("[ ] Performing seminario calculations")
    bonds, angles, dummy_dihedrals = multi_seminario(metal_charge, metal_name, keywords=keywords, mult=mult,
                                                     improper_metal=improper_metal, donors=donors)
    print("\t[ ] Seminario finished, now coping")
    for a, topol in enumerate(topols):
        topols[a] = copy_bonds(topols[a], bonds[a], metal_name)
        topols[a] = copy_angles(topols[a], angles[a], metal_name)
        topols[a] = copy_dihedrals(topols[a], dummy_dihedrals[a], metal_name)
    print("\t[+] Copied!")

    print("[ ] Calculating charges")
    partial_charges = calculate_charges(metal_charge, metal_name, vdw_data_name=vdw_type, mult=mult)

    print("\t[ ] Copying charges")
    for n_site, topol in enumerate(topols):
        for idx, atom in enumerate(topol.atoms):
            atom.charge = partial_charges[n_site][idx]
    print("\t[+] Charges calculated !")

    print("[ ] Cutting extra atoms")
    File = open("INFO.dat")  # TODO this here is wierd
    text = File.read()
    File.close()

    for line in text.splitlines():
        if "extra_atoms:" in line:
            extra_atoms = list(map(int, line[12:].split(',')))

    for a, topol in enumerate(topols):
        topol.write(f"old_new_topol{a:d}.top")

        topol.strip(f"@{','.join(list(map(str, np.array(extra_atoms) + 1))):s}")
        topol.write(f"new_topol{a:d}.top")
        new_cage = MDAnalysis.Universe(f"site{a:d}.pdb")
        new_cage.atoms.write(f"old_new_template{a:d}.pdb")

        print(f"We use this on site{a:d}.pdb:")
        print(f"not index {' '.join(list(map(str, np.array(extra_atoms)))):s}")
        new_cage.select_atoms(f"not index {' '.join(list(map(str, np.array(extra_atoms)))):s}").write(
            f"template{a:d}.pdb")

    print("Hydrogens removed!")

    print("Program finishess successfully!")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-o", help="Clean structure")
    parser.add_argument("-metal_name", help="Name of the metal")
    parser.add_argument("-metal_charge", help="Charge of the metal")
    parser.add_argument("-keywords", help="keywords for QM", nargs='+')
    parser.add_argument("-vdw_type", default='uff',
                        help="Type of parameters for VdW (available: uff, merz-tip3p, merz-opc3, merz-spc/e, merz-tip3p-fb, merz-opc, merz-tip4p-fb, merz-tip4-ew, zhang-tip3p, zhang-opc3, zhang-spc/e, zhang-spc/eb, zhang-tip3p-fb, zhang-opc, zhang-tip4p/2005, zhang-tip4p-d, zhang-tip4p-fb, zhang-tip4p-ew")
    parser.add_argument("-mult", default=1, help="multiplicity")
    parser.add_argument("-improper_metal", action='store_true', default=False,
                        help="Calculate improper dihedral of the metal-aromatic (default:False)")
    parser.add_argument("-donors", nargs='+', default=['N', 'S', 'O'],
                        help="Donors from the connected ligands, usually electronegative atom, such as N, S, O, but sometimes metal is connected to carbon", )
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    # print("AAAA", args.keywords)
    parametrize(args.f, args.metal_name, int(args.metal_charge), args.keywords, vdw_type=args.vdw_type,
                mult=int(args.mult), improper_metal=args.improper_metal, donors=args.donors)
