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
    adjust_charge, adjust_pair_exclusions

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
        print("directory", os.getcwd())
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

            mapping_fp_to_new, _ = find_mapping_of_fingerprint_on_metal_and_its_surroundings(cage_coord, site.index, site.metal_name, site.fp_syst, cutoff=site.ligand_cutoff)

            self.topol_new = adjust_charge(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_bonds(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_angles(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_dihedrals(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_impropers(self.topol_new, site.fp_topol, mapping_fp_to_new)

            self.topol_new = adjust_pair_exclusions(self.topol_new, site.fp_topol, mapping_fp_to_new)

        logger.info('Finished')

    ''' TODO remove (2023/07/04)
    def create_cage_and_linker_cgbind(self, smiles, arch_name, metal, metal_charge):
        
        # Creates cage and linker using cgbind. Then it paramterizes the linker using antechamber
        # 
        # :param smiles:
        # :param arch_name:
        # :param metal:
        # :param metal_charge:
        # :return:
        

        try:
            from cgbind import Linker, Cage
            from antechamber_interface import antechamber
        except:
            raise

        logger.info("[ ] Calling cgbind to create cage")
        linker = Linker(smiles=smiles, arch_name=arch_name)
        cage = Cage(linker, metal=metal, metal_charge=metal_charge)



        cage.print_xyz_file('cage.xyz')
        linker.print_xyz_file('linker.xyz')
        syst = MDAnalysis.Universe('linker.xyz')
        syst.atoms.write('linker.pdb')
        logger.info("[ ] Calling antechamber to parametrize linker") # TODO, that should not be hidden here
        antechamber('linker.pdb', 'linker.top')
    '''

    ''' TODO remove (2023/07/04)
    def load_cage(self, cage_file, ligand_file):
        
        Copies the cage file to the class. Renumbers the ligands to match the MD topology. Extracts some properties

        :param cage_file:
        :param ligand_file:
        :return:
        
        #cage = MDAnalysis.Universe(cage_file)
        #cage.atoms.write("temp.gro")
        #self.crystal2pdb("temp.gro", ligand_file, "cage.gro", metal_name=self.metal_name)

        #self.cage = MDAnalysis.Universe("cage.gro")
        #self.ligand = pmd.load_file(ligand_file)

        # Find metal indeces:
        
        # self.metal_indices = [a for a, name in enumerate(self.cage.atoms.names) if
        #                  name[:len(self.metal_name)].upper() == self.metal_name]
        # self.n_metals = len(self.metal_indices)
        


        logger.info(self.metal_indices)
        logger.info(self.metal_name)
        logger.info(f"[test {(len(self.cage.atoms) - self.n_metals) % len(self.ligand.atoms):} {len(self.cage.atoms):} "
                    f"{self.n_metals:} {len(self.ligand.atoms):}")
        assert (len(self.cage.atoms) - self.n_metals) % len(self.ligand.atoms) == 0
        self.n_ligands = int((len(self.cage.atoms) - self.n_metals) / len(self.ligand.atoms))


        # Find the lowest distance between metal sites, if there are more than 1
        if self.n_metals>1:
            mm_distances = distance_array(self.cage.atoms[self.metal_indices].positions, self.cage.atoms[self.metal_indices].positions)
            self.m_m_cutoff = np.min(mm_distances[mm_distances>0.1])-0.1 #this acctually does not make it better TODO
            #self.m_m_cutoff=9
        else:
            self.m_m_cutoff = 1e10


        logger.info(f"The distance between closest metals: {self.m_m_cutoff:f}")
        logger.info(f"Detected M{self.n_metals:d}L{self.n_ligands:d} cage")

        # Find how many ligands are bound to single metal:
        metal_index = self.metal_indices[0] # TODO check if all have the same number of ligands bound to single side
        '''

    ''' TODO remove (2023/07/04)
    def load_fingerprint_old(self, cage_file): #TODO remove
        
        # Loads files from the library into the class, if the name is not known, it will try to guess the fingerprint
        # 
        # :param name_of_binding_side:
        # :return:
        


        if self.name_of_binding_side is None:
            self.name_of_binding_side = guess_fingerprint(cage_file, self.metal_indices[0], metal_name=self.metal_name, fingerprint_guess_list=None,
                              m_m_cutoff=10)
        #else:
        #    self.name_of_binding_side = name_of_binding_side  # TODO check if it exists


        # Input
        self.fp_topol, self.fp_syst= load_fp_from_file(self.name_of_binding_side, self.fingerprint_style)

        metal_position = self.fp_syst.atoms[0].position
        nometal_position = self.fp_syst.atoms[1:].positions
        self.ligand_cutoff = np.min([np.max(distance_array(metal_position, nometal_position)) + 2.0, self.m_m_cutoff])


        #self.topol_fp = pmd.load_file(f'{os.path.dirname(__file__):s}/library/{self.name_of_binding_side:s}.top')
        #self.syst_fingerprint = MDAnalysis.Universe(f"{os.path.dirname(__file__):s}/library/{self.name_of_binding_side:s}.pdb")
        return True
    '''

    ''' TODO remove (2023/07/04) 
    def crystal2pdb(self, crystal_pdb, topology_itp, output, metal_name=""):
        
        # Function which renumbers crystal_pdb that it maches topology of ligand, it adds metals at the begining, creates
        #   the file output with renumbered atoms
        #
        # :param crystal_pdb:
        # :param topology_itp:
        # :param output:
        # :param metal_name:
        # :return:
        

        crystal = MDAnalysis.Universe(crystal_pdb)
        ligand = pmd.load_file(topology_itp)

        topology = MDAnalysis.Universe.empty(len(ligand.atoms), trajectory=True)
        topology.add_TopologyAttr('name')
        for a in range(len(topology.atoms)):
            topology.atoms[a].name = ligand.atoms[a].name

        table = MDAnalysis.topology.tables.vdwradii
        table['Cl'] = table['CL']

        for atom in crystal.atoms:
            atom.name = atom.name.upper()
        metal_name = metal_name.upper()

        Gtop = nx.Graph([(bond.atom1.idx, bond.atom2.idx) for bond in ligand.bonds])
        nx.set_node_attributes(Gtop, {atom.idx: atom.name[0] for atom in ligand.atoms}, "name")

        if len(metal_name) > 0:
            metal_ions = crystal.select_atoms("name {:s}*".format(metal_name))
            ligands = crystal.select_atoms("not name {:s}*".format(metal_name))
            G1 = nx.Graph(
                MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table,  box=crystal.dimensions))
            nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")
            new_ligands = [metal_ions]
        else:
            ligands = crystal.atoms
            G1 = nx.Graph(
                MDAnalysis.topology.guessers.guess_bonds(ligands.atoms, ligands.atoms.positions, vdwradii=table,  box=crystal.dimensions))
            nx.set_node_attributes(G1, {atom.index: atom.name[0] for atom in ligands.atoms}, "name")
            new_ligands = []

        for z, nodes in enumerate(nx.connected_components(G1)):
            new_ligand = topology.copy()
            #print(nodes)
            print(ligands)
            print(list(nodes))
            print(list(map(str,list(nodes))))
            ligands.select_atoms(f" index {' '.join(list(map(str,list(nodes)))):}").write(f"temp.pdb")

            # logger.info(len(nodes))
            Gsub = G1.subgraph(nodes)

            #Sometimes guesser will not guess right the bonds, decreasing fudge factor might help...
            if len(Gsub)==len(Gtop) and len(Gsub.edges())!=len(Gtop.edges()):
                fudge_factor = 0.55
                while fudge_factor>0 and len(Gsub.edges())!=len(Gtop.edges()):
                    Gsub = nx.Graph(MDAnalysis.topology.guessers.guess_bonds(crystal.atoms[Gsub.nodes()], crystal.atoms[Gsub.nodes()].positions, vdwradii=table,
                                                                             box=crystal.dimensions, fudge_factor=fudge_factor))
                    nx.set_node_attributes(Gsub, {atom.index: atom.name[0] for atom in crystal.atoms[Gsub.nodes()]}, "name")
                    fudge_factor-=0.01
                if fudge_factor<0.1:
                    print("error, cannot find the guess for bonds")

            iso = isomorphism.GraphMatcher(Gsub, Gtop, node_match=lambda n1, n2: n1['name'] == n2['name'])
            if iso.is_isomorphic():
                for node in list(nodes):
                    # print(node, iso.mapping[node])
                    # print(
                    # print(crystal.atoms[node].name, topology.atoms[iso.mapping[node]].name)
                    new_ligand.atoms[iso.mapping[node]].position = crystal.atoms[node].position

                new_ligands.append(new_ligand.atoms)
                # print(len(new_ligands),new_ligands )

            else:
                logger.info("The subgraphs are not isomorphic, are you sure these are the same molecules?")
                logger.info(f"Number of atoms: from coordinates: {len(Gsub):d}, from_topology: {len(Gtop):d}")
                logger.info(f"Number of edges: from coordinates: {len(Gsub.edges()):d}, from_topology: {len(Gtop.edges()):d}")

                print(os.getcwd())
                print(list(nodes))
                #crystal.atoms[list(nodes)].write("temp.pdb")
                #new_ligand.atoms[list(Gsub.nodes())].write("temp.pdb")

                print("The subgraphs are not isomorphic, are you sure these are the same molecules?")
                print(f"Number of atoms: from coordinates: {len(Gsub):d}, from_topology: {len(Gtop):d}")
                print(f"Number of edges: from coordinates: {len(Gsub.edges()):d}, from_topology: {len(Gtop.edges()):d}")

                print(os.getcwd())
                raise Error

        new_cage = MDAnalysis.Merge(*new_ligands)
        new_cage.dimensions=crystal.dimensions
        new_cage.atoms.write(output)

    '''
    #def prepare_new_topology(self, cage_coord, cage_topol, metal_name, metal_charge, ligand_file=None):
    def prepare_new_topology(self, cage_coord, cage_topol):
        '''
        metal = pmd.load_file(f'{os.path.dirname(__file__):s}/library/M.itp')
        metal.atoms[0].name = self.metal_name
        metal.atoms[0].charge = self.metal_charge
        metal.atoms[0].mass = name2mass[self.metal_name.title()]

        cage_topol = metal * self.n_metals + self.ligand * self.n_ligands
        cage_topol.write('test_cage.top', [list(range(self.n_ligands + self.n_metals))])
        '''


        #self.n_ligands, self.metal_indices, self.n_metals = prepare_initial_topology(cage_file, metal_name, metal_charge, "cage.gro", "test_cage.top", ligand_topol=ligand_file) # TODO can this output be nicer?

        self.cage = MDAnalysis.Universe(cage_coord)
        #logger.info(f"[ ] Created initial topol with {self.n_ligands:} ligands, {self.n_metals} metals (indecies: {self.metal_indices:})")

        self.topol_new = pmd.load_file(cage_topol, parametrize=False)  # , skip_bonds=True)
        # we need to copy paramters back, for some reasons if paramtrize is False
        # parmed is sometimes wierd :-(

        topol_new2 = pmd.load_file(cage_topol)
        for a, _ in enumerate(self.topol_new.atoms):
            self.topol_new.atoms[a].type = topol_new2.atoms[a].type
            self.topol_new.atoms[a].epsilon = topol_new2.atoms[a].epsilon
            self.topol_new.atoms[a].sigma = topol_new2.atoms[a].sigma
            self.topol_new.atoms[a].rmin = topol_new2.atoms[a].rmin
            self.topol_new.atoms[a].charge = topol_new2.atoms[a].charge

        # Find the lowest distance between metal sites, if there are more than 1

        ''' # TODO this should moe
        if self.n_metals > 1:
            mm_distances = distance_array(self.cage.atoms[self.metal_indices].positions, self.cage.atoms[self.metal_indices].positions)
            self.m_m_cutoff = np.min(mm_distances[mm_distances>0.1])-0.1 #this acctually does not make it better TODO
        else:
            self.m_m_cutoff = 1e10
            #self.m_m_cutoff=9
        '''
    ''' TODO remove (2023/07/04)
    def adjust_charge(self, mapping_fp_to_new):
        logger.info("   [ ] Changing charges and atomtypes")
        sum_of_charge_diffrences = 0

        for a in mapping_fp_to_new:
            logger.info(f"          {self.topol_new.atoms[mapping_fp_to_new[a]].type:s} "
                        f"{self.topol_new.atoms[mapping_fp_to_new[a]].name:s} "
                        f"{self.topol_new.atoms[mapping_fp_to_new[a]].charge:} --> "
                        f"{self.fp_topol.atoms[a].type, self.fp_topol.atoms[a].name:} {self.topol_new.atoms[mapping_fp_to_new[a]].charge + self.fp_topol.atoms[a].charge:}")


            atom = self.fp_topol.atoms[a]

            if atom.type not in self.topol_new.parameterset.atom_types.keys():
                logger.info(f"      [^] Adding new atomtype: {atom.type:s}")
                atomtype = pmd.topologyobjects.AtomType(name=atom.type, number=atom.number, mass=atom.mass)
                self.topol_new.parameterset.atom_types[atom.type] = atomtype

            self.topol_new.atoms[mapping_fp_to_new[a]].type = self.fp_topol.atoms[a].type
            self.topol_new.atoms[mapping_fp_to_new[a]].epsilon = self.fp_topol.atoms[a].epsilon
            self.topol_new.atoms[mapping_fp_to_new[a]].sigma = self.fp_topol.atoms[a].sigma
            self.topol_new.atoms[mapping_fp_to_new[a]].rmin = self.fp_topol.atoms[a].rmin
            self.topol_new.atoms[mapping_fp_to_new[a]].charge += self.fp_topol.atoms[a].charge  # topol_fp.atoms[a].charge

            sum_of_charge_diffrences += self.fp_topol.atoms[a].charge

    def adjust_bonds(self, mapping_fp_to_new):
        logger.info("   [ ] Adding new bonds to topology")
        for bond_fp in self.fp_topol.bonds:
            found = False
            for bond_new in self.topol_new.bonds:

                # Check if atoms in mapping (that they are not the additional atoms)
                if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                    if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                         mapping_fp_to_new[bond_fp.atom2.idx]) or
                            ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                              mapping_fp_to_new[bond_fp.atom1.idx]))):
                        if (bond_fp.type != bond_new.type):
                            logger.info(f"      [o] Diffrent bond type  {bond_new.type:} -> {bond_fp.type:}")
                            logger.info(f"          Mapping ({bond_fp.atom1.idx},{bond_fp.atom2.idx}) to ({bond_new.atom1.idx:},{bond_new.atom2.idx:})")
                            logger.info(f"          Mapping ({bond_fp.atom1.name},{bond_fp.atom2.name}) to ({bond_new.atom1.name:},{bond_new.atom2.name:})")
                            type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                                          list=self.topol_new.bond_types)
                            #if type_to_assign not in self.topol_new.bond_types:
                            self.topol_new.bond_types.append(type_to_assign)

                            #bond_new.type = deepcopy(bond_fp.type)
                            bond_new.type = type_to_assign #bond_fp.type

        logger.info("   [ ] Adding new bonds to topology")

        for bond_fp in self.fp_topol.bonds:
            found = False

            for bond_new in self.topol_new.bonds:
                # if bond_fp.atom1.name=="ZN":
                #    print(bond_fp.atom1.idx , bond_new)

                # print(bond_new, bond_new.atom1.name)
                # Check if atoms in mapping (that they are not the additional atoms)
                if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                    # print("It is in the mapping")

                    if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                         mapping_fp_to_new[bond_fp.atom2.idx]) or
                            ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                              mapping_fp_to_new[bond_fp.atom1.idx]))):
                        if (bond_fp.type != bond_new.type):
                            logger.info(f"      [o] Diffrent bond type {bond_new.type:} {bond_fp.type:}")
                        found = True
                else:

                    found = True # TODO not sure why this is here

            if not found:
                # type_to_assign= bond_fp.type
                # type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq)
                # print("AAA", bond_fp.type.k, bond_fp.type.req, len(topol_new.bond_types))

                type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                              list=self.topol_new.bond_types)

                logger.info(f"      [o] New bond: {mapping_fp_to_new[bond_fp.atom1.idx]:} "
                            f"{mapping_fp_to_new[bond_fp.atom2.idx]:}, {type_to_assign:}")
                # if type_to_assign not in topol_new.bond_types:
                self.topol_new.bond_types.append(type_to_assign)

                atom1 = self.topol_new.atoms[mapping_fp_to_new[bond_fp.atom1.idx]]
                atom2 = self.topol_new.atoms[mapping_fp_to_new[bond_fp.atom2.idx]]
                self.topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))


    def adjust_angles(self, mapping_fp_to_new):
        logger.info("   [ ] Adding new angles to topology")
        for angle_fp in self.fp_topol.angles:
            found = False
            if angle_fp.atom1.idx in mapping_fp_to_new and angle_fp.atom2.idx in mapping_fp_to_new and angle_fp.atom3.idx in mapping_fp_to_new:
                for angle_new in self.topol_new.angles:
                    #Check if atoms in mapping (that they are not the additional atoms)
                    if ((angle_new.atom1.idx==mapping_fp_to_new[angle_fp.atom1.idx] and angle_new.atom2.idx==mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx==mapping_fp_to_new[angle_fp.atom3.idx]) or
                        ((angle_new.atom1.idx==mapping_fp_to_new[angle_fp.atom3.idx] and angle_new.atom2.idx==mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx==mapping_fp_to_new[angle_fp.atom1.idx]))):
                        if (angle_fp.type!=angle_new.type):
                            logger.info("      [o] Diffrent angle type {angle_new.type:} -> {angle_fp.type:}")
                            #angle_new.funct = angle_fp.funct
                            #angle_new.type.k = angle_fp.type.k
                            #angle_new.type.theteq = angle_fp.type.theteq


                            type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq, list = self.topol_new.angle_types)

                            #if type_to_assign not in topol_new.angle_types:
                            self.topol_new.angle_types.append(type_to_assign)
                            angle_new.type = type_to_assign

                        found = True
            else:
                    logger.info(f"[-] Not found in fingerprint: {angle_fp.atom1.idx+1:d} {angle_fp.atom2.idx+1:d} {angle_fp.atom3.idx+1:d}")
                    logger.info(f"              {self.fp_topol.atoms[angle_fp.atom1.idx]:}")
                    logger.info(f"              {self.fp_topol.atoms[angle_fp.atom2.idx]:}")
                    logger.info(f"              {self.fp_topol.atoms[angle_fp.atom3.idx]:}")
                    logger.info(f"But it is standard bond so it should be there")
                    logger.info(f"And you should be worry if it is not the end of fingerprint")
                    found = True
                    #raise

            if not found:
                #type_to_assign= angle_fp.type
                type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq,  list = self.topol_new.angle_types )
                #if type_to_assign not in topol_new.angle_types:
                self.topol_new.angle_types.append(type_to_assign)
                logger.info(f"new type {len(self.topol_new.angle_types):}")

                logger.info(f"      [x] New angle:{angle_fp.atom1.name:}-{angle_fp.atom2.name:}-{angle_fp.atom3.name:} "
                            f"{mapping_fp_to_new[angle_fp.atom1.idx]:} {mapping_fp_to_new[angle_fp.atom2.idx]:} "
                            f"{mapping_fp_to_new[angle_fp.atom3.idx]:} {type_to_assign:}")
                #print("              ", topol_fp.atoms[angle_fp.atom1.idx])
                #print("              ", topol_fp.atoms[angle_fp.atom2.idx])
                #print("              ", topol_fp.atoms[angle_fp.atom3.idx])
                atom1= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom1.idx]]
                atom2= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom2.idx]]
                atom3= self.topol_new.atoms[mapping_fp_to_new[angle_fp.atom3.idx]]

                self.topol_new.angles.append(pmd.topologyobjects.Angle(atom1, atom2, atom3, type=type_to_assign))

    def adjust_dihedrals(self, mapping_fp_to_new):
        logger.info("   [ ] Adding new dihedrals to topology")
        for dihedral_fp in self.fp_topol.dihedrals:
            found = False
            if dihedral_fp.atom1.idx in mapping_fp_to_new and dihedral_fp.atom2.idx in mapping_fp_to_new and dihedral_fp.atom3.idx in mapping_fp_to_new and dihedral_fp.atom4.idx in mapping_fp_to_new:
                for dihedral_new in self.topol_new.dihedrals:
                    #Check if atoms in mapping (that they are not the additional atoms)
                    if (((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom1.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom4.idx]) or
                        ((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom4.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom1.idx]))) and
                       (dihedral_new.type.per==dihedral_fp.type.per)):
                        if (dihedral_fp.type!=dihedral_new.type):
                            logger.info(f"      [a] Diffrent Dihedral type: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:}")
                            logger.info(f"           {dihedral_new.type:}")
                            logger.info(f"        -> {dihedral_fp.type:}")

                            type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedral_fp.type.phi_k, per=dihedral_fp.type.per,
                                                                             phase=dihedral_fp.type.phase, scee=dihedral_fp.type.scee,
                                                                             scnb=dihedral_fp.type.scnb, list = self.topol_new.dihedral_types)
                            #if type_to_assign not in topol_new.dihedral_types:
                            self.topol_new.dihedral_types.append(type_to_assign)
                            dihedral_new.type = type_to_assign

                        found = True
            else:
                    logger.info(f"[-] Not found in finger print: {dihedral_fp.atom1.idx+1:} {dihedral_fp.atom2.idx+1:} {dihedral_fp.atom3.idx+1:} {dihedral_fp.atom4.idx+1:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom1.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom2.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom3.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom4.idx]:}")
                    logger.info("But it is standard so it should be there")
                    logger.info("And you should be worry if it is not the end of fingerprint")
                    found = True # TODO I don't remember why this is True
                    #raise
                    
                    # for dihedral in self.topol_fp.dihedrals:
                    # if dihedral.atom1.idx==0:
                    # print(dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1)
                    

            if not found:
                type_to_assign= pmd.topologyobjects.DihedralType(phi_k=dihedral_fp.type.phi_k, per=dihedral_fp.type.per,
                                                                phase=dihedral_fp.type.phase, scee=dihedral_fp.type.scee,
                                                                scnb=dihedral_fp.type.scnb, list = self.topol_new.dihedral_types)
                #if type_to_assign not in topol_new.dihedral_types:
                self.topol_new.dihedral_types.append(type_to_assign)
                logger.info(f"      [b] New dihedral: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:} "
                            f"{mapping_fp_to_new[dihedral_fp.atom1.idx]:} {mapping_fp_to_new[dihedral_fp.atom2.idx]:} "
                            f"{mapping_fp_to_new[dihedral_fp.atom3.idx]:} {type_to_assign:}")

                atom1 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom1.idx]]
                atom2 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom2.idx]]
                atom3 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom3.idx]]
                atom4 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom4.idx]]

                self.topol_new.dihedrals.append(pmd.topologyobjects.Dihedral(atom1, atom2, atom3, atom4, type=type_to_assign))


    def adjust_impropers(self, mapping_fp_to_new):
        # This is almost the same as above, with the diffrence on being improper dihedrals

        logger.info("   [ ] Adding new improper dihedrals to topology")
        for dihedral_fp in self.fp_topol.impropers:
            found = False
            if dihedral_fp.atom1.idx in mapping_fp_to_new and dihedral_fp.atom2.idx in mapping_fp_to_new and dihedral_fp.atom3.idx in mapping_fp_to_new and dihedral_fp.atom4.idx in mapping_fp_to_new:
                for dihedral_new in self.topol_new.impropers:
                    #Check if atoms in mapping (that they are not the additional atoms)
                    if (((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom1.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom4.idx]) or
                        ((dihedral_new.atom1.idx==mapping_fp_to_new[dihedral_fp.atom4.idx] and dihedral_new.atom2.idx==mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom3.idx==mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom4.idx==mapping_fp_to_new[dihedral_fp.atom1.idx]))) and
                       (dihedral_new.type.per==dihedral_fp.type.per)):
                        if (dihedral_fp.type!=dihedral_new.type):
                            logger.info(f"      [a] Diffrent Dihedral type: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:}")
                            logger.info(f"           {dihedral_new.type:}")
                            logger.info(f"        -> {dihedral_fp.type:}")

                            type_to_assign = pmd.topologyobjects.ImproperType(psi_k=dihedral_fp.type.psi_k, psi_eq=dihedral_fp.type.psi_eq,
                                                                              list = self.topol_new.improper_types)
                            #if type_to_assign not in topol_new.dihedral_types:
                            self.topol_new.improper_types.append(type_to_assign)
                            dihedral_new.type = type_to_assign

                        found = True
            else:
                    logger.info(f"[-] Not found in finger print: {dihedral_fp.atom1.idx+1:} {dihedral_fp.atom2.idx+1:} {dihedral_fp.atom3.idx+1:} {dihedral_fp.atom4.idx+1:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom1.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom2.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom3.idx]:}")
                    logger.info(f"               {self.fp_topol.atoms[dihedral_fp.atom4.idx]:}")
                    logger.info("But it is standard so it should be there")
                    logger.info("And you should be worry if it is not the end of fingerprint")
                    found = True # TODO I don't remember why this is True
                    #raise
                    
                    # for dihedral in self.topol_fp.dihedrals:
                    # if dihedral.atom1.idx==0:
                    # print(dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1)
                    

            if not found:
                type_to_assign = pmd.topologyobjects.ImproperType(psi_k=dihedral_fp.type.psi_k,
                                                                  psi_eq=dihedral_fp.type.psi_eq,
                                                                  list=self.topol_new.improper_types)

                #if type_to_assign not in topol_new.dihedral_types:
                self.topol_new.improper_types.append(type_to_assign)
                logger.info(f"      [b] New dihedral: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:} "
                            f"{mapping_fp_to_new[dihedral_fp.atom1.idx]:} {mapping_fp_to_new[dihedral_fp.atom2.idx]:} "
                            f"{mapping_fp_to_new[dihedral_fp.atom3.idx]:} {type_to_assign:}")

                atom1 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom1.idx]]
                atom2 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom2.idx]]
                atom3 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom3.idx]]
                atom4 = self.topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom4.idx]]

                self.topol_new.impropers.append(pmd.topologyobjects.Improper(atom1, atom2, atom3, atom4, type=type_to_assign))
    '''

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
