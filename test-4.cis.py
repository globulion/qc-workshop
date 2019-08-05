#!/usr/bin/python3
import tutor.psithon.util
from tutor.project_3.cis import CIS

mol = tutor.psithon.util.psi_molecule_from_file('water.xyz')

set {
 scf_type df
 basis 6-31G*
 e_convergence 1e-8
 reference uhf
}

cis = CIS.create(mol, verbose=True)
cis.run()
