Tutorials
=====

.. _tutorials:


Tutorial 1: Force-field parameters from SMILES string
------------

It is straightforward to obtain force-field parameters just from SMILES string by combining the functionality of metallicious with `cgbind <https://github.com/duartegroup/cgbind/tree/master>`_:

.. code-block:: python

    from cgbind import Linker, Cage
    from metallicious import supramolecular_structure

    # Create coordination file from cgbind
    linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')
    cage = Cage(linker, metal='Pd')
    cage.print_xyz_file(filename="cage.xyz")

    # Parametrize the cage, use AmberTools to parametrize organic linkers
    cage = supramolecular_structure('cage.xyz', metal_charges={'Pd': 2}, vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

Tutorial 2: Parametrization of template for gallium cages (known issue and solution)
------------

Metallicious uses educated guesses along the way. While we observe that this works in most cases, we also found one outlier: parametrization of the template for [Ga4L6]12- cages.
To parametrize the template, *metalliciuos* uses the charge of the metal (user input) and guessed charge of the attached ligand fragment.
The charge of the organic ligand is guessed by rdkit from the 3D structure of the atoms. While this works for most of the molecules,
sometimes there might be ambiguity about the structure. In this example, the coordinating ligand can have two chemical structures:

.. image:: images/gallium.png
  :width: 200
  :align: center
  :alt:

*metallicious* incorrectly guesses the neutral charge of the ligand. However, this behaviour can be corrected by overwriting guessed charges:

.. code-block:: python

    from metallicious import supramolecular_structure
    cage = supramolecular_structure('nonbonded.pdb', {'Ga': (3, 1)},
                                    topol = 'noncovalent_complex.top',
                                    vdw_type='uff', search_library=False)
    cage.extract_unique_metal_sites()
    cage.unique_sites[0].ligand_charges = [-2, -2, -2] #this specifies ligands' charges
    cage.parametrize()

While it is worth noting this, this template is part of the default library, so it is already parametrized.


Tutorial 3: Using custom template
------------

One might decide that they parametrized the template using a specific technique and did not save it to the library. One might use this template on purpose for specific application

.. code-block:: python

    cage = supramolecular_structure(f'cage.xyz',
                                    metal_charges={'Pd': 2},
                                    vdw_type = 'uff',
                                    search_library=False)
    cage.extract_unique_metal_sites() # extracts template structures (needed to detect sites)
    cage.sites[0].fp_coord_file = f'selected_template.pdb' # loads coordinates of the template
    cage.sites[0].fp_topol_file = f'selected_template.top' # loads force-field parameters of the template
    cage.sites[0].load_fingerprint() # loads template into the class
    cage.sites[0].set_cutoff() # reads the radius of the template
    cage.unique_sites = [] # overwrites list with templates to parameterize (no template parametrization needed)
    cage.parametrize(out_coord='saturated_template.pdb', out_topol='saturated_template.top')


Tutorial 4: Custom template library
------------

As default *metallicious* has a library of the templates parametrized using D3BJ-PBE0/def2-SVP and they are saved in directory of *metallicious* in library folder.
One might however opt for higher level of theory or include implicit solvent effect. The easiest is to create new directory for library of templates.
This requires only change of the library_directory, as there are not templates inside, we need to parametrize them using new method:


.. code-block:: python

    from metallicious import supramolecular_structure
    cage = supramolecular_structure('nonbonded.pdb', {'Pd': (2, 1)},
                                    topol = 'noncovalent_complex.top',
                                    vdw_type='uff'
                                    library_path=f'/path/to/new/library/',
                                    keywords = ['CPCM(Water)', 'PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq'])
    cage.parametrize()