#!/usr/bin/python3
"""
A Psi4 input script to compute CIS energy from a SCF reference

Algorithms were developed for ilustrative purposes.
"""

import time
import numpy
import psi4
from abc import ABC, abstractmethod
from ..psithon.util import two_index_transform, four_index_transform


class CIS(ABC):
  """
 Implementation of CIS method.
"""
  reference_types = ['rhf', 'uhf']

  def __init__(self, mol, verbose, save_states):
      ABC.__init__(self)
      self._common_init(mol, verbose, save_states)


  # ---> public interface <--- #

  @classmethod
  def create(cls, mol, verbose=True, save_states=None, reference='rhf'):#TODO
      raise NotImplementedError

  def run(self):
      "Run CIS calculations"
      self._run_scf()
      self._prepare_for_cis()
      self._build_hamiltonian()
      self._diagonalize()


  # ---> protected interface <--- #

  def _common_init(self, mol, verbose, save_states):
      self.mol = mol
      self.verbose = verbose
      self.save_states = save_states
      #
      self.ref_wfn = None
      self.scf_e = None
      self.e_0 = None
      self.nuclear_repulsion_energy = mol.nuclear_repulsion_energy()
      self.N = self.nuclear_repulsion_energy
      #
      self.hamiltonian = None
      self.E = None
      self.W = None
      #
      self.nmo = None
      self.naocc = None
      self.nbocc = None
      self.navir = None
      self.nbvir = None
      self.ndet = None
      #
      self.Ca_occ = None
      self.Cb_occ = None
      self.Ca_vir = None
      self.Cb_vir = None
      #
      self.eri_OVOV = None
      self.eri_OOVV = None
      self.eri_OVov = None
      self.eri_oovv = None
      self.eri_ovov = None
      #
      self.same_ab = None
      self.ci_l = None

  def _run_scf(self):#TODO
      raise NotImplementedError

  def _prepare_for_cis(self):#TODO
      raise NotImplementedError

  @abstractmethod
  def _set_beta(self, H, eri): pass

  @abstractmethod
  def _set_scf_reference(self): pass

  def _build_hamiltonian(self):#TODO
      raise NotImplementedError

  def _diagonalize(self):#TODO
      raise NotImplementedError
    
     

class RCIS(CIS):
  def __init__(self, mol, verbose, save_states):
      CIS.__init__(self, mol, verbose, save_states)
      self.same_ab = True

  def _set_scf_reference(self):
      psi4.core.set_global_option('reference', 'rhf')

  def _set_beta(self, H, eri):#TODO
      raise NotImplementedError


class UCIS(CIS):
  def __init__(self, mol, verbose, save_states):
      CIS.__init__(self, mol, verbose, save_states)
      self.same_ab = False

  def _set_scf_reference(self):
      psi4.core.set_global_option('reference', 'uhf')

  def _set_beta(self, H, eri):#TODO
      raise NotImplementedError
