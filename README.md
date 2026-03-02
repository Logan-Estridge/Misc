# Miscellaneous scripts for my research

## fig_all_kde_using_class.py
A tool for generating Kernel Density Estimate (KDE) plots for PLUMED CV data.

## sort_atoms_paba.py
Given an xyz file of pABA with atoms tagged as belonging to "clusters" (from Ovito, export with cluster info in the first column), this script consistently orders the atoms within each molecule and separates all the clusters. Useful for PLUMED analysis where having separate molecules with a consistent atom ordering is crucial.
