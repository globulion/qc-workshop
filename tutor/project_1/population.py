#!/usr/bin/python3
"""
 Module for population analyses
"""
import psi4
import numpy
from ..psithon.util import matrix_power

__all__ = ["atomic_charges", "Loc"]

def atomic_charges(wfn, kappa=0.0):
    """
 Compute atomic partial charges as a function of kappa parameter. 
 
 Input:
   wfn   - psi4.core.Wavefunction object
   kappa - parameter in the interval [0, 1]

 Returns:
   numpy.ndarray of shape (wfn.molecule().natom(), ) with partial charges [A.U.]
 
 Notes:
   o kappa = 0 corresponds to Mulliken charges
   o kappa = 1/2 corresponds to Lowdin charges
"""
    assert(0.0<=kappa<=1.0), "Kappa must be between 0.0 and 1.0"
    mol = wfn.molecule()
    bfs = wfn.basisset()
    # initialize partial charges with atomic numbers
    charges = numpy.array([mol.Z(x) for x in range(mol.natom())],numpy.float64)
    # compute bond order matrix
    P = wfn.Da().clone()
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

class Loc:
  def __init__(self, wfn, method='BOYS'):
      self.wfn = wfn
      self.method = method
      localizer_a = psi4.core.Localizer.build(method, wfn.basisset(), wfn.Ca_subset("AO","OCC"))
      localizer_b = psi4.core.Localizer.build(method, wfn.basisset(), wfn.Cb_subset("AO","OCC"))
      localizer_a.localize()
      localizer_b.localize()
      self.La = localizer_a.L
      self.Lb = localizer_b.L
      self.Ua = localizer_a.U
      self.Ub = localizer_b.U
  def lmoc(self):
      "Compute LMO centroids" 
      # initialize
      lmoc_a = numpy.zeros((self.wfn.nalpha(), 3), dtype=numpy.float64)
      lmoc_b = numpy.zeros((self.wfn.nbeta (), 3), dtype=numpy.float64)
      # compute dipole integrals
      mints = psi4.core.MintsHelper(self.wfn.basisset())
      D = mints.ao_dipole()
      # compute centroids
      for z in range(3):
          la_z = -numpy.einsum("ai,bi,ab->i", self.La, self.La, D[z])
          lb_z = -numpy.einsum("ai,bi,ab->i", self.Lb, self.Lb, D[z])
          lmoc_a[:,z] = la_z
          lmoc_b[:,z] = lb_z
      return lmoc_a, lmoc_b
  def el_charges(self):
      qa = -numpy.ones(self.wfn.nalpha(), dtype=numpy.float64)
      qb = -numpy.ones(self.wfn.nbeta (), dtype=numpy.float64)
      return qa, qb
  def el_dipoles(self): 
      da = 0; db=0
      da = -numpy.zeros((self.wfn.nalpha(),3), dtype=numpy.float64)
      db = -numpy.zeros((self.wfn.nbeta (),3), dtype=numpy.float64)
      return da, db
  def el_quadrupoles(self):
      Qa = -numpy.zeros((self.wfn.nalpha(),6), dtype=numpy.float64)
      Qb = -numpy.zeros((self.wfn.nbeta (),6), dtype=numpy.float64)
      # compute quadrupole integrals
      mints = psi4.core.MintsHelper(self.wfn.basisset())
      D = mints.ao_quadrupole()
      # compute quadrupoles
      for z in range(6):
          la_z = +numpy.einsum("ai,bi,ab->i", self.La, self.La, D[z])
          lb_z = +numpy.einsum("ai,bi,ab->i", self.Lb, self.Lb, D[z])
          Qa[:,z] = la_z
          Qb[:,z] = lb_z
      la, lb = self.lmoc()
      for i in range(self.wfn.nalpha()):
          ri2 = numpy.outer(la[i], la[i]) 
          ri2xx = ri2[0,0]
          ri2xy = ri2[0,1]
          ri2xz = ri2[0,2]
          ri2yy = ri2[1,1]
          ri2yz = ri2[1,2]
          ri2zz = ri2[2,2]
          r2 = numpy.array([ri2xx,ri2xy,ri2xz,ri2yy,ri2yz,ri2zz])
          Qa[i] += r2
      for i in range(self.wfn.nbeta()):
          ri2 = numpy.outer(lb[i], lb[i]) 
          ri2xx = ri2[0,0]
          ri2xy = ri2[0,1]
          ri2xz = ri2[0,2]
          ri2yy = ri2[1,1]
          ri2yz = ri2[1,2]
          ri2zz = ri2[2,2]
          r2 = numpy.array([ri2xx,ri2xy,ri2xz,ri2yy,ri2yz,ri2zz])
          Qb[i] += r2
      return Qa, Qb
