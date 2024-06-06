Examples
=====

.. _examples:

The examples here are

Example 1
------------

Simple paramterization of the cage with 2 metals. As an imput we provide coordination file in PBD format and non-bonded
force-field paramters in *.top format (GROMACS). Morover LJ paramters of the Ru and Pd are taken from universal force-field
(vdw_type=uff), as they are only one available for Ru metal.

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.pdb', metal_charges={'Ru': 2, 'Pd':2 }, topol='ru_pd.top', vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')


Example 2
------------

In case of lack of the topology simple parametrization using AmberTools with GAFF force-field is done of the linker.
*metallicious* will prepare the non-bonded force-field paramters of the whole structure which is saved to folder "init_topol".

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

The parameters can be saved to AMBER format force-field

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, vdw_type='uff')
    cage.parametrize(out_coord='out.inpcrd', out_topol='out.prmtop', prepare_initial_topology=True)


Example 3
------------

If the cage is homoleptic and you have available force-field paramters of the linker, you can use them as an input to
generate initial topology of the whole structure.

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, vdw_type='uff')
    cage.prepare_initial_topology(homoleptic_ligand_topol='linker.top')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')

Introduction to examples 4 and 5:
------------
Cage shown below have rather complex metal binding site.

The template for this metal site is not available in metallicious. Therfore running command:

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

raises an error "Topology file not specified, please provide topology, or use prepare_initial_topology=True".

Two solutions are available:
1. parametrize the metal site. This is done by providing multiplicity of the metal:
2. use truncation scheme (caution needed)

Example 4: Parametrization of new template
------------

If paramters for the template are not available, you might decide to parametrize them. In *metallicious* this is done by
specifing the multiplicity of the metal which also enembales QM calculation. For this functionality, the additional
dependencies (see installation guide) are needed (`autode <https://github.com/duartegroup/autodE>`_, `ORCA <https://orcaforum.kofo.mpg.de/app.php/portal>`_, and`psiRESP <https://github.com/lilyminium/psiresp>`_).

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

Bare in mind that parametrization of template is time-consuming. It will perform DFT optimalisation using ORCA/autodE. By default autodE uses 4 cores but this can be modified:

.. code-block:: python
    from metallicious import supramolecular_structure
    import autode as ade
    ade.Config.n_cores = 8
    cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

By default parametrization is done on D3BJ-PBE0/def2-SVP (keywords = ['PBE0', 'D3BJ', 'def2-SVP', 'tightOPT', 'freq']).
This can be changed by specifing keywords in class supramolecular structure class:

.. code-block:: python
    from metallicious import supramolecular_structure
    cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc', keywords= ['B3LYP', '6-31G*', 'tightOPT', 'freq'])
    cage.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)


Example 5
------------

Truncation scheme are very simple

.. code-block:: python
    from metallicious import supramolecular_structure
    # This will not work becasue there is no exact template for this site:
    # cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc')

    cage = supramolecular_structure('cage.pdb', topol='topol.top', metal_charges={'Pd':2 }, vdw_type='merz-opc', truncation_scheme='dihedral')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')

