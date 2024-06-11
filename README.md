<p align="center">
<img src="images/logo.png" width="100"/>
</p>

# Metallicious 

Metallicious (a playful combination of "metal" and "delicious") is an automated tool for creating force fields for 
metal-containing systems with a covalent model of the metal. By utilizing a library of templates, metallicious identifies 
the template that matches the metal site in the structure. It copies the bonded parameters from the template and performs 
charge redistribution to account for charge transfer. In cases where no suitable template is found, metallicious 
automatically performs parameterization.

<img src="images/summary.png" height="150"/>

Metallicious works with minimal user input, relying heavily on educated guesses, which may not always yield the expected 
results. Therefore, it is recommended to use the tool with caution.



## Installation:
```
conda install rdkit autode psiresp mdanalysis networkx qcelemental==0.25.1 ambertools --channel conda-forge
pip install metallicious
```
## Quick start
Parametrization of structure with coordinates saved as `supramolecular_cage.xyz` (*.xyz, *.pdb, *gro, etc. formats 
supported by MDAnalysis) with (nonbonded) force-field parameters `supramolecular_cage.top` (*top, *prmtop, etc. supported by ParmEd), 
which consists of two metals Pd2+ and Ru2+ and organic linkers, the input file might look like this:

```
from metallicious import supramolecular_structure
cage = supramolecular_structure('supramolecular_cage.xyz',
                                metal_charges={'Ru': 2, 'Pd':2 },
                                topol='supramolecular_cage.top',
                                vdw_type='uff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
```
The `supramolecular_structure` function takes a coordination file, topology file , dictionary of metal ions along with 
their charge (and in case template parametrization is needed multiplicity) and type of Lennard-Jones library as an input.

See [examples/](https://github.com/tkpiskorz/metallicious/tree/main/metallicious/examples) for more examples and 
[metallicious.readthedocs.io](https://metallicious.readthedocs.io/en/latest/) for additional documentation.

## Command line
It is also possible to use the metallicious just form command line. For example:
```
metallicious -f supramolecular_cage.xyz -vdw_type merz-tip3p -metal_and_charges Pd 2 Ru 2 -prepare_topol
```
For details, see:
```
metallicious -h
```
