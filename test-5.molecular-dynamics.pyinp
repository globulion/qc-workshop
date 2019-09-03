#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_3.cis import CIS
import tutor.project_4.surface_hopping as sh, numpy, time
from tutor.project_1.population import atomic_charges, Loc
numpy.set_printoptions(precision=4, linewidth=200, suppress=True)

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

set {
 scf_type df
 basis sto-3g
 e_convergence 1e-8
 reference rhf
 guess read
 puream 1
}

#system = sh.DynamicalSystem(mol, dt_class=0.2, qm_method='psi4scf', nstates=1, init_state=0, seed=404, temperature=300.0)
system = sh.DynamicalSystem(mol, dt_class=0.2, qm_method='mycis', nstates=8, init_state=3, seed=4, temperature=0.0)
system.run(100, center_mode='qm')
