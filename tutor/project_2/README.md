# Hartree-Fock Model

The usual starting point in the description of the electronic ground state
of molecules is through the Hartree-Fock, or molecular orbital, approximation.
Here we consider HF theory for the closed shell case,
where the total energy of the molecule is given by

<img src="../../doc/figures/equations/energy.png" height="40"/>

where the effective one-electron Fock matrix **F** is defined 

<img src="../../doc/figures/equations/fock-matrix.png" height="40"/>

**D**, the one-particle density matrix,

<img src="../../doc/figures/equations/density-matrix.png" height="40"/>

is build up from the occupied block of the molecular orbitals LCAO-MO matrix
**C**, given by the generalized eigenvalue problem

<img src="../../doc/figures/equations/orbital-energies.png" height="40"/>

# Implementation


-----------
[Main Page](https://github.com/globulion/qc-workshop)
