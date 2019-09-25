#!/usr/bin/python3
"""Quantum Computer Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from .ciwavefunction import HF_CIWavefunction, CIS_CIWavefunction
from ..project_2.cis import CIS
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
  def _compute_energy(self):#TODO
  def _compute_forces(self):#TODO

class ExcitedState_Computer(Computer):
  def __init__(self, molecule, nstates):
      Computer.__init__(self, molecule)
      self.nstates = nstates

  def _compute_forces(self):#TODO


class CIS_Computer(ExcitedState_Computer):
  def __init__(self, molecule, nstates):
      ExcitedState_Computer.__init__(self, molecule, nstates) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule, nstates):
      CIS_Computer.__init__(self, molecule, nstates)

  def _compute_energy(self):#TODO
