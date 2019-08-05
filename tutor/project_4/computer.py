#!/usr/bin/python3
"""Quantum Computer Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from .ciwavefunction import HF_CIWavefunction, CIS_CIWavefunction
from ..project_3.cis import CIS
from ..psithon.util import rearrange_eigenpairs, check_sim

__all__ = ["Computer"]

class Computer(ABC):
  def __init__(self, molecule):
      ABC.__init__(self)
      self.restart_wfn = None
      #
      self.molecule = molecule
      self.nuclear_repulsion_energy = molecule.nuclear_repulsion_energy()
      #
      self.ciwfn = None
      self.nstates = None
      self.forces = None

  @classmethod
  def create(cls, method, molecule, nstates=1):
      m = method.lower()
      if m == "psi4scf": return psi4SCF_Computer(molecule)
      elif m == "mycis": return myCIS_Computer(molecule, nstates)
      else: raise ValueError("Wrong method chosen for computer")

  def update(self, xyz):
      self.molecule.set_geometry(xyz)
      self.nuclear_repulsion_energy = self.molecule.nuclear_repulsion_energy()

  def compute(self): 
      self.ciwfn = self._compute_energy()
      self.restart_wfn = self.ciwfn.ref_wfn
      self.forces = self._compute_forces()

  @abstractmethod
  def _compute_energy(self): pass

  @abstractmethod
  def _compute_forces(self): pass



class psi4SCF_Computer(Computer):
  def __init__(self, molecule): 
      Computer.__init__(self, molecule)
      self.nstates = 1
  def _compute_energy(self):
      #if self.restart_wfn is not None:
      #   e, w = psi4.energy("SCF", molecule=self.molecule, ref_wfn=self.restart_wfn, return_wfn=True)
      #else:
      #   e, w = psi4.energy("SCF", molecule=self.molecule, return_wfn=True)
      e, w = psi4.energy("SCF", molecule=self.molecule, return_wfn=True)
      E = numpy.array([e - self.molecule.nuclear_repulsion_energy()])
      state = HF_CIWavefunction(w, E)
      return state
  def _compute_forces(self):
      g = -psi4.gradient("SCF", molecule=self.molecule).to_array(dense=True)
      return numpy.array([g])

class ExcitedState_Computer(Computer):
  def __init__(self, molecule, nstates):
      Computer.__init__(self, molecule)
      self.nstates = nstates

  def _compute_forces(self):
      xyz = self.molecule.geometry()
      d = 0.0001
      e_0 = self.ciwfn.ci_e + self.nuclear_repulsion_energy
      w_0 = self.ciwfn.ci_c.copy()
      n_states = len(e_0)

      g = numpy.zeros((n_states, self.molecule.natom(), 3), numpy.float64)
      for a in range(self.molecule.natom()):
          for x in range(3):                                                           
              xyz_xa1 = xyz.clone()                                  
              xyz_xa1.set(a, x, xyz.get(a, x) + d)
              #
              self.molecule.set_geometry(xyz_xa1)
              s = self._compute_energy()
              e = s.ci_e + self.molecule.nuclear_repulsion_energy()
              w = s.ci_c
              e, w, sim = rearrange_eigenpairs(w, w_0, e, return_sim=True)
              log = check_sim(sim)
              if 'ERROR' in log: raise ValueError("Error in numerical differentiation of forces! %s" % sim)
              g_xa1 = (e - e_0) / d
              #
              g[:,a,x] = g_xa1
      self.molecule.set_geometry(xyz)
      #f = -numpy.array(f).reshape(self.molecule.natom(),3,n_states).transpose(2,0,1)
      return -g


class CIS_Computer(ExcitedState_Computer):
  def __init__(self, molecule, nstates):
      ExcitedState_Computer.__init__(self, molecule, nstates) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule, nstates):
      CIS_Computer.__init__(self, molecule, nstates)
  def _compute_energy(self):
      cis = CIS.create(self.molecule, verbose=False, save_states=self.nstates-1, reference='rhf')
      cis.run()
      state = CIS_CIWavefunction(cis.ref_wfn, cis.E, cis.W)
      return state
