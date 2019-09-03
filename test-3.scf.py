#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_2.scf import SCF
import numpy

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

set {
 scf_type pk
 basis 6-31G*
 e_convergence 1e-11
 reference rhf
 puream True
}

# Psi4 result
e_hf, w_hf = energy('scf', molecule=mol, return_wfn=True)
print(" Energy Psi4 result = %14.6f [a.u.]" % e_hf)

# Guess: H-core
print("\n Running SCF for Hcore guess")
scf_1 = SCF(mol)
scf_1.run(verbose=True, guess=None)

# Guess: scaled Fock matrix
print("\n Running SCF for scaled Fock guess")
scf_2 = SCF(mol)
scf_2.run(verbose=True, guess=w_hf.Fa().to_array()*0.2)
