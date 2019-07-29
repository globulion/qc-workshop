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
  def _overlap_between_slater_determinants(self, other):
      # overlap between AO's
      mints = psi4.core.MintsHelper(self.bfs)
      s_ao  = mints.ao_overlap(self.bfs, other.bfs) 
      #
      S_nm  = numpy.zeros((self.ndet, other.ndet), numpy.float64)
      # sum over all Slater determinants
      for n in range(self.ndet):
          cai, cbi = self._construct_c(n)
          for m in range(other.ndet):
              caj, cbj = other._construct_c(m)
              D_ij_a = numpy.linalg.multi_dot([cai.T, s_ao, caj])
              D_ij_b = numpy.linalg.multi_dot([cbi.T, s_ao, cbj])
              S_nm[n,m] = numpy.linalg.det(D_ij_a) * numpy.linalg.det(D_ij_b)
      return numpy.abs(S_nm)
  def _construct_c(self, n):
      det = self.ci_l[n]
      ca  = self.ca_o.copy()
      cb  = self.cb_o.copy()
      if det.is_reference: return ca, cb
      if det.is_single:
         if det.change_alpha: ca[:,det.rule[0]] = self.ca_v[:,det.rule[1]]
         else:                cb[:,det.rule[0]] = self.cb_v[:,det.rule[1]]
      if det.is_double: raise NotImplementedError
      return ca, cb

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
  def make_ci_l(self):
      dets = []
      dets.append(Reference_SlaterDeterminant(self.naocc, self.nbocc, self.nmo))
      for i in range(self.naocc):
          for a in range(self.navir):
              rule = (i,a)
              det = Single_SlaterDeterminant(self.naocc, self.nbocc, self.nmo, rule)
              dets.append(det)
      for i in range(self.nbocc):
          for a in range(self.nbvir):
              rule = (-i,-a)
              det = Single_SlaterDeterminant(self.naocc, self.nbocc, self.nmo, rule)
              dets.append(det)
      return dets



