
Example of usage:

From smiles:
```
cgbind2pmd -smiles "C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1" -arch_name m2l4 -metal pd -charge 2
```

From cage coordination file and topology of the linker:
```commandline
cgbind2pmd -f cage.pdb -linker_topol gromacs.top -metal Pd -metal_charge 2
```




