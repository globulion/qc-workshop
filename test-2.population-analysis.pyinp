#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_1.population import atomic_charges, Loc
import numpy
numpy.set_printoptions(precision=4, linewidth=200, suppress=True)

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

set {
 scf_type df
 basis 6-31++G(d,p)
 e_convergence 1e-8
}

e, w = energy('scf', molecule=mol, return_wfn=True)

# atomic charges from SCF density
for kappa in numpy.linspace(0.0, 0.5, 20):
    q = atomic_charges(w, kappa)
    print(" kappa = %5.2f   q(O) = %14.4f   q(H) = %14.4f" % (kappa, q[0], q[1]))

# LMO's and distributed quadrupole moments
loc=  Loc(w)
print(loc)
