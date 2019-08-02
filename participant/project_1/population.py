#!/usr/bin/python3
"""
 Module for population analyses
"""
import psi4
import numpy
from ..psithon.util import matrix_power

__all__ = ["atomic_charges", "Loc"]

def atomic_charges(wfn, kappa=0.0):#TODO
    raise NotImplementedError

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
  def lmoc(self):#TODO
      "Compute LMO centroids" 
      # return alpha and beta LMO centroids as (la, lb) where la is a numpy.ndarray of shape (nmo, 3)
      raise NotImplementedError

  def el_charges(self):
      qa = -numpy.ones(self.wfn.nalpha(), dtype=numpy.float64)
      qb = -numpy.ones(self.wfn.nbeta (), dtype=numpy.float64)
      return qa, qb

  def el_dipoles(self): 
      da = -numpy.zeros((self.wfn.nalpha(),3), dtype=numpy.float64)
      db = -numpy.zeros((self.wfn.nbeta (),3), dtype=numpy.float64)
      return da, db

  def el_quadrupoles(self):#TODO
      # return distributed alpha and beta quadrupole moments as (Qa, Qb) where Qa is numpy.ndarray of shape (nmo, 6)
      # Hint: there are only 6 independent Cartesian components of quadrupole moment.
      #       In Psi4, the turn of indices is XX, XY, XZ, YY, YZ, ZZ
      #                                        0   1   2   3   4   5
      raise NotImplementedError

  def __repr__(self):
      "Print LMO centroids and LMO quadrupole moments"
      la, lb = self.lmoc()
      Qa, Qb = self.el_quadrupoles()
      # convert units to angs and debye*angs
      la*= psi4.constants.bohr2angstroms
      lb*= psi4.constants.bohr2angstroms
      Qa*= psi4.constants.bohr2angstroms*psi4.constants.dipmom_au2debye
      Qb*= psi4.constants.bohr2angstroms*psi4.constants.dipmom_au2debye
      #
      log = " \n"
      log+= " Localized Orbital Centroids and Quadrupole Moments (Loc = %s)\n" % self.method.upper()
      log+= " Orbital  r_i [Angs]                         Q_i [Debye Angs]\n" 
      log+= " <alpha>\n"
      for i in range(self.wfn.nalpha()):
          log += "%2da   " % (i+1)
          log += " %8.3f %8.3f %8.3f     " % tuple(la[i])
          log += " %11.3f %11.3f %11.3f  " % tuple(Qa[i,:3])
          log += "  <aver>=%8.3f\n" % ((Qa[i,0]+Qa[i,3]+Qa[i,5])/3.0)
          log += " "*38
          log += " %11.3f %11.3f %11.3f   \n" % (Qa[i,1], Qa[i,3], Qa[i,4])
          log += " "*38
          log += " %11.3f %11.3f %11.3f   \n" % (Qa[i,2], Qa[i,4], Qa[i,5])
          log += "\n"
      log+= " <beta>\n"
      for i in range(self.wfn.nbeta ()):
          log += "%2db   " % (i+1)
          log += " %8.3f %8.3f %8.3f     " % tuple(lb[i])
          log += " %11.3f %11.3f %11.3f  " % tuple(Qb[i,:3])
          log += "  <aver>=%8.3f\n" % ((Qb[i,0]+Qb[i,3]+Qb[i,5])/3.0)
          log += " "*38
          log += " %11.3f %11.3f %11.3f   \n" % (Qb[i,1], Qb[i,3], Qb[i,4])
          log += " "*38
          log += " %11.3f %11.3f %11.3f   \n" % (Qb[i,2], Qb[i,4], Qb[i,5])
          log += "\n"

      return str(log)
