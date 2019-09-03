#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_2.scf import SCF
import numpy, psi4

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

psi4.set_options({"scf_type"     : "pk", 
                  "basis"        : "6-31G*",
                  "e_convergence": 1e-11,
                  "puream"       : True})
psi4.core.set_output_file('test-3.out', True)

# Psi4 result
e_hf, w_hf = psi4.energy('scf', molecule=mol, return_wfn=True)
print(" Energy Psi4 result = %14.6f [a.u.]" % e_hf)

# Guess: H-core
print("\n Running SCF for Hcore guess")
scf_1 = SCF(mol)
scf_1.run(verbose=True, guess=None)

# Guess: scaled Fock matrix
print("\n Running SCF for scaled Fock guess")
scf_2 = SCF(mol)
scf_2.run(verbose=True, guess=w_hf.Fa().to_array()*0.2)
