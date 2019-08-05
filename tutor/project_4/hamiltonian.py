#!/usr/bin/python3
"""Quantum Computer Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import psi4
import numpy
from .computer import Computer

__all__ = ["Isolated_Hamiltonian"]

class Hamiltonian(ABC):
  def __init__(self, method, aggregate, nstates, method_low=None):
      ABC.__init__(self)
      self.aggregate = aggregate
      self.computers = []
      self.computers.append(Computer.create(method, aggregate.qm, nstates))
      for i in range(1,aggregate.nfrags-1):
          self.computers.append(Computer.create(method_low, aggregate.bath[i]))

  @classmethod
  def create(cls, molecule, method, method_low=None): pass

  def compute(self):
      "Solve Schrodinger Equation for Extended Molecular Aggregate"
      # unperturbed Hamiltonian
      self.computers[0].compute()
      # environment
      if self.aggregate.nfrags > 1:
         self.iterate_bath() 
      return
 
  @abstractmethod
  def iterate_bath(self): pass

class Isolated_Hamiltonian(Hamiltonian):
  def __init__(self, molecule, method, nstates):
      Hamiltonian.__init__(self, molecule, method, nstates, method_low=None)
  def iterate_bath(self): raise ValueError("This Hamiltonian is isolated")

class QMMM_Hamiltonian(Hamiltonian):
  "MM potential is set on L"

class QMQM_Hamiltonian(Hamiltonian):
  "QM potential is set on L"

class EE_QMQM_Hamiltonian(QMQM_Hamiltonian):
  "Instantaneous point charges from wavefunction on L are QM potential"

class PE_QMQM_Hamiltonian(QMQM_Hamiltonian):
  "Polarization embedding potential is set on L"

class PDE_QMQM_Hamiltonian(PE_QMQM_Hamiltonian):
  "Polarization density embedding potential is set on L"
