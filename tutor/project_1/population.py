#!/usr/bin/python3
"""
 Module for population analyses
"""
import numpy
from ..psithon.util import matrix_power

def atomic_charges(wfn, kappa=0.0):
    """
 Compute atomic partial charges as a function of kappa parameter. 
 
 Input:
   wfn   - psi4.core.Wavefunction object
   kappa - parameter in the interval [0, 1/2]

 Returns:
   numpy.ndarray of shape (wfn.molecule().natom(), ) with partial charges [A.U.]
 
 Notes:
   o kappa = 0 corresponds to Mulliken charges
   o kappa = 1/2 corresponds to Lowdin charges
"""
    assert(0.0<=kappa<=0.5), "Kappa must be between 0.0 and 0.5"
    mol = wfn.molecule()
    bfs = wfn.basisset()
    # initialize partial charges with atomic numbers
    charges = numpy.array([mol.Z(x) for x in range(mol.natom())],numpy.float64)
    # compute bond order matrix
    P = wfn.Da()
    P.add(wfn.Db())
    # extract one-electron integrals
    S = wfn.S()
    # compute SkPSk1 matrix
    Sk = matrix_power(S, kappa)
    Sk1= matrix_power(S, 1.0 - kappa)
    G = numpy.linalg.multi_dot([Sk, P, Sk1])
    # 
    for i in range(bfs.nbf()):
        x = bfs.function_to_center(i)
        charges[x] -= G[i,i]
    return charges
