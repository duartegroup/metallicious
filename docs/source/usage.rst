Usage
=====

.. _usage:

Python module
------------

Quick_start

Parametrization of structure with coordinates saved as `supramolecular_cage.xyz` with (nonbonded) topology `supramolecular_cage.top` (of the whole structure):


    from metallicious import supramolecular_structure
    cage = supramolecular_structure('supramolecular_cage.xyz',
                                    metal_charges={'metal name 1': charge of metal 1(integer), 'metal name 2':charge of metal 2(integer),...},
                                    topol='supramolecular_cage.top', vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')


For example, for the structure ru_pd.xyz with force-field parameters saved as ru_pd.top, which consists of two metals Pd2+ and Ru2+, the input file looks like this:


    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, topol='ru_pd.top', vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')



## Initial topology file
If you don't have a topology file, you can generate a simple force-field parametrization using General Amber Force-field (GAFF):

    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, vdw_type='uff')
    cage.prepare_initial_topology()
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')

However, we do not intend to automate the parametrization of the organic part of the molecule.
Please refer to specialized tools such as gromacs, atb, ambertools, and charmm-gui.


# Handling missing templates

The number of combination of possible ligands and metal results that inevitibly you will encounter metal site for which there is no template.
In such cases to solutions are possible:

### 1. Parametrize template
We recommend to run template parametrization on HPC/cluster as it can take some time (our experience is ~8h on 8 CPUs per template).

Specifying explicitly the metal multiplicity using the metal_charge_mult variable instead of metal_charges, will automatically inform metallicious to be ready to parametrize the template


    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charge_mult = {'Ru': (2,1), 'Pd':(2,1)}, vdw_type='uff')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')


### 2. Truncate existing template

If an exact template is unavailable in the library, you can truncate part of an existing template.
Truncation is based on the distance from the metal centre, such as 4-bonds away ("dihedral"), 3-bonds away ("angles"), or 1-bond away ("bonds").
Such a strategy is fast but results in a loss of accuracy.

For example:

    from metallicious import supramolecular_structure
    cage = supramolecular_structure('ru_pd.xyz', metal_charge_mult = {'Ru': (2,1), 'Pd':(2,1)}, truncation_scheme = 'dihedral')
    cage.parametrize(out_coord='out.pdb', out_topol='out.top')

The extended list of supramolecular_structure class:

    :param filename: (str) name of the coordination file
    :param metal_charge_mult:  (dict) the names charges, and multiplicity of the metals in format
                                    {metal_name: (metal_charges, multiplicity)}
    :param metal_charges: (dict) the names and charges of metals in the input structure in format:
                                      {metal_name1: metal_charges1, metal1_name2: metal_charge2}
    :param vdw_type: (str) name of LJ dataset used for metal paramters
    :param topol: (str) path to topology (optional)
    :param keywords: list(str) the keywords for QM calculations
    :param improper_metal: (bool) if True it will parametrize the improper dihedral involving metal
    :param donors: (list(str)) list of atom elements with which metal forms bond
    :param library_path: (str) directory of template library, be default where the script is
    :param ff: (str) parametrization protocol for small organic molecules (only gaff available)
    :param search_library: (bool) if True, metallicious searrched templates in template library,
                if False, it will parametrize template
    :param fingerprint_guess_list: (list(str)) list of templates to check
    :param truncation_scheme: (str) name of the truncation scheme
    :param covalent_cutoff: (float) if metal-atoms smaller then cutoff it creates bonds ligand with metal


Bash command line
------------

It is also possible to use the metallicious just form command line. For example:

metallicious -f cage.xyz -vdw_type merz-tip3p -metal_and_charges Pd 2 -prepare_topol

For details, see:

metallicious -h

.. list-table:: Title
   :widths: 25 25 50
   :header-rows: 1

   * - [copy from SI]
     - Heading row 1, column 2
     - Heading row 1, column 3
   * - Row 1, column 1
     -
     - Row 1, column 3
   * - Row 2, column 1
     - Row 2, column 2
     - Row 2, column 3


-f Metaloorganic structure (*.gro, *.pdb, etc. all supported by MDAnalysis)
-p Metaloorganic topology of the whole structure (*.top, *.prmtop, etc. all supported by ParmEd)", default=False)
-of", help="Output metaloorganic structure (*.gro, *.pdb, etc. supported by MDAnalysis)", default='out.pdb')
parser.add_argument("-op", help="Output metaloorganic topology (*.top, *.prmtop, etc. supported by ParmEd)", default='out.top')

parser.add_argument("-metal_and_charge",nargs='+', help="Metal names and charges (optionally, multiplicity, when parametrization needed), for example: Pd 2 1 Ru 2 1")
parser.add_argument("-keywords", help="autodE keywords for QM calculations (default: PBE0 D3BJ def2-SVP tightOPT freq)", nargs='+')
parser.add_argument("-LJ_type", default='merz-opc',
                    help="Type of parameters for Lennard-Jones paramters (default: merz-opc; available: uff, merz-tip3p, merz-opc3, merz-spc/e, merz-tip3p-fb, merz-opc, merz-tip4p-fb, merz-tip4-ew, zhang-tip3p, zhang-opc3, zhang-spc/e, zhang-spc/eb, zhang-tip3p-fb, zhang-opc, zhang-tip4p/2005, zhang-tip4p-d, zhang-tip4p-fb, zhang-tip4p-ew")
parser.add_argument("-truncate", help="Truncation scheme (defualt: None; available: None, 3-bond (dihedral), 2-bond (angle), 1-bond(bond))", default=None)
parser.add_argument("-improper_metal", action='store_true', default=False,
                    help="Calculate improper dihedral of the metal-aromatic (default: False)")
parser.add_argument("-donors", nargs='+', default=['N', 'S', 'O'],
                    help = "Donors from the connected ligands, usually electronegative atom, such as N, S, O, but sometimes metal is connected to carbon (default: N S O)")
parser.add_argument("-prepare_initial_topol", action='store_true', default=False, help="Prepare initial topology using GAFF")
parser.add_argument("-linker_topol", default=None, help="Linker force-field (topology) parameters")
return parser.parse_args()



Available parameters
------------

Default templates

Lennard-Jones