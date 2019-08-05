#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_2.scf import SCF
import numpy

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

set {
 scf_type df
 basis 6-31G*
 e_convergence 1e-8
 reference rhf
}

scf_1 = SCF(mol)
scf_1.run(verbose=True, guess=None)

e_hf, w_hf = energy('scf', molecule=mol, return_wfn=True)
print(" Energy Psi4 result = %14.6f [a.u.]" % e_hf)
