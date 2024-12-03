# Flexible Protein-Ligand Docking Paper's files

This repository contains input simulations files from 
["Interactive molecular dynamics in virtual  reality for accurate flexible protein-ligand docking"](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0228461).
which are:

1. ** 2g5n_complex.xml ** : trypsin and indole-amidine, based on PDB file 2G5N. 
    An 800 kJ/mol/nm2 restraint was added to all backbone atoms of trypsin.
2.  ** hiv1_complex.xml ** : HIV-1 protease and amprenavir, based on PDB file 1HPV.
    An 800 kJ/mol/nm2 was applied to all backbone atoms of HIV-1 complex, excluding those which 
    make up the flaps which gate the active site (defined as residues 49 to 55 in chain A and 
    residues 48 to 54 in chain B).
3.  ** 6hcx_complex.xml ** : neuraminidase and zanamivir, based on PDB file 6HCX.
    No restraint is applied to the complex.

All files have been converted from NarupaXR (.xml) to Nanover iMD-VR (.xml)