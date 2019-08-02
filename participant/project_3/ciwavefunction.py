#!/usr/bin/python3
"""CI Wavefunction Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import psi4
import numpy

__all__ = ["SlaterDeterminant",
              "Reference_SlaterDeterminant",
              "Single_SlaterDeterminant",
           "CIWavefunction",
              "HF_CIWavefunction",
              "CIS_CIWavefunction",]

class SlaterDeterminant(ABC):
  def __init__(self, nao, nbo, nmo, rule):
      ABC.__init__(self)
      self.is_reference = False
      self.is_single    = False
      self.is_double    = False
      self.is_triple    = False
      self.rule = rule
      self.nao = nao
      self.nbo = nbo
      self.nav = nmo - nao
      self.nbv = nmo - nbo
      self.nmo = nmo

class Reference_SlaterDeterminant(SlaterDeterminant):
  def __init__(self, nao, nbo, nmo):
      SlaterDeterminant.__init__(self, nao, nbo, nmo, rule=())
      self.is_reference = True

class Single_SlaterDeterminant(SlaterDeterminant):
  def __init__(self, nao, nbo, nmo, rule):
      SlaterDeterminant.__init__(self, nao, nbo, nmo, rule)
      self.is_single = True
      self.change_alpha = True if self.rule[0] > 0 else False


class CIWavefunction(ABC):
  "Multiconfigurational state composed of Slater determinants as basis set"
  def __init__(self, ref_wfn, E, W):
      ABC.__init__(self)
      self.ref_wfn = ref_wfn
      self.ci_e = E.copy()
      self.ci_c = W.copy()
      self.ca_o = ref_wfn.Ca_subset("AO","OCC").to_array(dense=True)
      self.cb_o = ref_wfn.Cb_subset("AO","OCC").to_array(dense=True)
      self.ca_v = ref_wfn.Ca_subset("AO","VIR").to_array(dense=True)
      self.cb_v = ref_wfn.Cb_subset("AO","VIR").to_array(dense=True)
      self.bfs  = ref_wfn.basisset()
      self.naocc= ref_wfn.nalpha()
      self.nbocc= ref_wfn.nbeta()
      self.nmo  = ref_wfn.nmo()
      self.navir= self.nmo - self.naocc
      self.nbvir= self.nmo - self.nbocc
      self.ci_l = self.make_ci_l()
      self.ndet = len(self.ci_l)

  @abstractmethod
  def make_ci_l(self): pass

  def overlap(self, other):
      assert(self.ndet == other.ndet)
      s_nm  = self._overlap_between_slater_determinants(other)
      S_IJ  = numpy.linalg.multi_dot([self.ci_c.T, s_nm, other.ci_c])
      return S_IJ

  def _overlap_between_slater_determinants(self, other):#TODO
      raise NotImplementedError

  def _construct_c(self, n):#TODO
      raise NotImplementedError

class HF_CIWavefunction(CIWavefunction):
  def __init__(self, ref_wfn, E):
      W = numpy.array([[1.0]])
      CIWavefunction.__init__(self, ref_wfn, E, W)
  def make_ci_l(self):
      dets = []
      dets.append(Reference_SlaterDeterminant(self.naocc, self.nbocc, self.nmo))
      return dets


class CIS_CIWavefunction(CIWavefunction):
  def __init__(self, ref_wfn, E, W):
      CIWavefunction.__init__(self, ref_wfn, E, W)
  def make_ci_l(self):#TODO
      raise NotImplementedError

