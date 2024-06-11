<p align="center">
<img src="images/logo.png" width="100"/>
</p>

# Metallicious 

Metallicious (a playful combination of "metal" and "delicious") is an automated tool for creating force fields for metal-containing systems with a covalent model of the metal. By utilizing a library of templates, Metallicious identifies the template that matches the metal site in the structure. It copies the bonded parameters from the template and performs charge redistribution to account for charge transfer. In cases where no suitable template is found, Metallicious automatically performs parameterization.

<img src="images/summary.png" height="150"/>

Metallicious works with minimal user input, relying heavily on educated guesses, which may not always yield the expected results. Therefore, it is recommended to use the tool with caution.



## Installation:
```
conda install rdkit autode psiresp mdanalysis networkx qcelemental==0.25.1 ambertools --channel conda-forge
pip install metallicious
```
## Quick start
Parametrization of structure with coordinates saved as `supramolecular_cage.xyz` with (nonbonded) topology `supramolecular_cage.top` (of the whole structure): 
```
from metallicious import supramolecular_structure
cage = supramolecular_structure('supramolecular_cage.xyz',
                                metal_charges={'metal name 1': charge of metal 1(integer), 'metal name 2':charge of metal 2(integer),...},
                                topol='supramolecular_cage.top', vdw_type='uff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
```

For example, for the structure ru_pd.xyz with force-field parameters saved as ru_pd.top, which consists of two metals Pd2+ and Ru2+, the input file looks like this:
```
from metallicious import supramolecular_structure
cage = supramolecular_structure('ru_pd.xyz', metal_charges={'Ru': 2, 'Pd':2 }, topol='ru_pd.top', vdw_type='uff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top')
```
The `supramolecular_structure` function takes a coordination file (*.xyz, *.pdb, *.gro, etc. supported by MDAnalysis), topology file (*top, *prmtop, supported by ParmEd), dictionary of metal ions along with their charge (and in case template parametrization is needed multiplicity) and type of Van der Waals metal parameters as input.

See [examples/](https://github.com/tkpiskorz/metallicious/tree/main/metallicious/examples) for
more examples and [metallicious.readthedocs.io](https://metallicious.readthedocs.io/en/latest/) for
additional documentation.

## Command line
It is also possible to use the metallicious just form command line. For example:
```
metallicious -f cage.xyz -vdw_type merz-tip3p -metal_and_charges Pd 2 -prepare_topol
```
For details, see:
```
metallicious -h
```
