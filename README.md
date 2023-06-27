The CMP (covalent metal parametrizer) is automated too to derviation of metal force-field parameters for classical Molecular Dynamics and Molecular modeling. The cmp is based on finding patterns of metal binding sites in existing library. If paramters are not available in the library cmp invokes automatic procedure for automatic derivation of paramters.

Currently the code supports only organometalic structures with metals sepparated at least by 2 non metal atoms, that is structure cannot have metal clusters. Morover, all the same type metals have to have the same charge and multiplicity.

The code works on the minimial input from user.*  

Requirements:
* autode
* ORCA
* rdkit
* networkx
* MDAnalysis
* AmberTools

Installation:

Simple example:

```
from parametrize_new_sites import supramolecular_structure
cage = supramolecular_structure('ru_pd.xyz', {'Ru': (2,1), 'Pd':(2,1)}, vdw_type='uff')
cage.parametrize(out_coord='out.pdb', out_topol='out.top)
```

supramolecular_structure(filename, metal_dictionery) takes as input a coordination file (*.xyz, *pdb, *gro etc.) and dictionery of metal ions togheter with their charge and multiplicity.



*(Warning! Since the code has minimal input, there are many steps which might no go as expected. Please use with cautious.)

