#!/usr/bin/python3
"""
A Psi4 input script to compute CIS energy from a SCF reference

Algorithms were developed for ilustrative purposes.
"""

__authors__   = "Bartosz Błasiak"
__copyright__ = "(c) 2019"
__license__   = "BSD-3-Clause"
__date__      = "2019-07-23"

import time
import numpy
import psi4
from abc import ABC, abstractmethod
from ..psithon.util import two_index_transform, four_index_transform

#numpy.set_printoptions(precision=3, linewidth=200, suppress=True)

__all__ = [\
           "CIS",
              "RCIS",
              "UCIS",]



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
  def create(cls, mol, verbose=True, save_states=None, reference='rhf'):
      "Create CIS instance"
      if reference.lower() not in cls.reference_types:
         raise ValueError("Incorrect reference wavefunction type chosen. Only RHF and UHF are available")
      assert(not (reference.lower()=='rhf' and mol.multiplicity() != 1)), "RHF reference cannot be set for open-shell system!"
      # UCIS
      if mol.multiplicity()!=1 or reference.lower()=='uhf': return UCIS(mol, verbose, save_states)
      else: return RCIS(mol, verbose, save_states)

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
      #self.Da = None
      #self.Db = None
      #self.jk = None
      #
      self.same_ab = None
      self.ci_l = None

  def _run_scf(self):
      self._set_scf_reference()
      scf_e, wfn = psi4.energy('HF', molecule=self.mol, return_wfn=True)
      self.scf_e = scf_e
      self.e_0   = scf_e - self.nuclear_repulsion_energy
      self.ref_wfn = wfn

  def _prepare_for_cis(self):
      self.Ca_occ = self.ref_wfn.Ca_subset("AO","OCC")
      self.Ca_vir = self.ref_wfn.Ca_subset("AO","VIR")
      #self.Da = self.ref_wfn.Da()
      self.nmo = self.ref_wfn.nmo()
      self.naocc = self.ref_wfn.nalpha()
      self.navir = self.nmo - self.naocc
      #
      #self.jk = psi4.core.JK.build(self.ref_wfn.basisset(), jk_type='direct')
      #self.jk.set_memory(int(5e8))
      #self.jk.initialize()
      #
      mints = psi4.core.MintsHelper(self.ref_wfn.basisset())
      eri = numpy.asarray(mints.ao_eri())
      H = self.ref_wfn.H().to_array(dense=True)
      #
      Oa = self.Ca_occ.to_array(dense=True)
      Va = self.Ca_vir.to_array(dense=True)
      #
      #self.eri_OVOV = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Va, Oa, Va, eri)
      #self.eri_OOVV = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Oa, Va, Va, eri)
      self.eri_OVOV = four_index_transform(eri, Oa, Va, Oa, Va)
      self.eri_OOVV = four_index_transform(eri, Oa, Oa, Va, Va)
      #
      #self.jk.C_clear()
      #self.jk.C_left_add(psi4.core.Matrix.from_array(self.Da, ""))
      #I = numpy.identity(self.Da.shape[0], numpy.float64)
      #self.jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      #self.jk.compute()
      #Ja= self.jk.J()[0].to_array(dense=True)
      #Ka= self.jk.K()[0].to_array(dense=True)
      ##
      #Ga = H+2.0*Ja-Ka
      Ga = self.ref_wfn.Fa().to_array(dense=True)
      #self.Fa_occ = numpy.einsum("ai,ab,bj->ij", Oa, Ga, Oa)
      #self.Fa_vir = numpy.einsum("ai,ab,bj->ij", Va, Ga, Va)
      self.Fa_occ = two_index_transform(Ga, Oa, Oa)
      self.Fa_vir = two_index_transform(Ga, Va, Va)
      #
      self._set_beta(H, eri)

  @abstractmethod
  def _set_beta(self, H, eri): pass

  @abstractmethod
  def _set_scf_reference(self): pass

  def _build_hamiltonian(self):
      # OO block
      self.hamiltonian[0, 0] = self.e_0
      # OS and SO blocks are zero
      None
      # SS block
      off_a = self.naocc*self.navir
      off_b = self.nbocc*self.nbvir
      for i in range(self.naocc):
          for a in range(self.navir):
              ia = self.navir*i + a
              # block AA
              for j in range(self.naocc):
                  for b in range(self.navir):
                      jb = self.navir*j + b
                      v = 0.0
                      if (i==j) and (a==b): v+= self.e_0
                      if (i==j): v+= self.Fa_vir[a,b]
                      if (a==b): v-= self.Fa_occ[i,j]
                      v += self.eri_OVOV[i,a,j,b] - self.eri_OOVV[i,j,a,b] 
                      #
                      self.hamiltonian[1+ia,1+jb] = v
              # block AB and BA
              for j in range(self.nbocc):
                  for b in range(self.nbvir):
                      jb = self.nbvir*j + b
                      v  = self.eri_OVov[i,a,j,b]
                      self.hamiltonian[1+ia,1+jb+off_a] = v
                      self.hamiltonian[1+jb+off_a,1+ia] = v
      if not self.same_ab:
         for i in range(self.nbocc):                                                         
             for a in range(self.nbvir):
                 ia = self.nbvir*i + a
                 # block BB
                 for j in range(self.nbocc):
                     for b in range(self.nbvir):
                         jb = self.nbvir*j + b
                         v = 0.0
                         if (i==j) and (a==b): v+= self.e_0
                         if (i==j): v+= self.Fb_vir[a,b]
                         if (a==b): v-= self.Fb_occ[i,j]
                         v += self.eri_ovov[i,a,j,b] - self.eri_oovv[i,j,a,b] 
                         #
                         self.hamiltonian[1+ia+off_a,1+jb+off_a] = v
      else:
          self.hamiltonian[(1+off_a):,(1+off_a):] = self.hamiltonian[1:1+off_a,1:1+off_a]
      del self.eri_ovov, self.eri_OVOV, self.eri_OVov, self.eri_OOVV, self.eri_oovv

  def _diagonalize(self):
      t = time.time()
      E, W = numpy.linalg.eigh(self.hamiltonian)
      if self.save_states is not None:
         E = E[  :(self.save_states+1)]
         W = W[:,:(self.save_states+1)]
      self.E = E
      self.W = W
      # 
      e_cis_ground = E[0] + self.nuclear_repulsion_energy
      if self.verbose: print('..finished diagonalization in %.3f seconds.\n' % (time.time() - t))
      if self.verbose: print('No of Determinants: % 16d' % (self.ndet))
      if self.verbose: print('SCF energy:         % 16.10f' % (self.scf_e))
      if self.verbose: print('CIS ground:         % 16.10f' % (e_cis_ground))
     
      hartree2eV = 27.211
      
      if self.verbose: print('\nCIS Excitation Energies:')
      if self.verbose: print('          Hartree                  eV')
      if self.verbose: print('--  --------------------  --------------------')
      for i in range(1, len(E)):
          excit_e = E[i] + self.nuclear_repulsion_energy - e_cis_ground
          if self.verbose: print('%2d %20.10f %20.10f' % (i, excit_e, excit_e * hartree2eV))

     
     

class RCIS(CIS):
  def __init__(self, mol, verbose, save_states):
      CIS.__init__(self, mol, verbose, save_states)
      self.same_ab = True

  def _set_scf_reference(self):
      psi4.core.set_global_option('reference', 'rhf')

  def _set_beta(self, H, eri):
      self.Cb_occ = self.Ca_occ
      self.Cb_vir = self.Ca_vir
      #self.Db = self.Da
      self.nbocc = self.naocc
      self.nbvir = self.navir
      self.ndet = 1 + self.naocc * self.navir + self.nbocc * self.nbvir
      self.hamiltonian = numpy.zeros((self.ndet, self.ndet),numpy.float64)
      #
      self.eri_ovov = self.eri_OVOV
      self.eri_OVov = self.eri_OVOV
      self.eri_oovv = self.eri_OOVV
      self.Fb_occ = self.Fa_occ
      self.Fb_vir = self.Fa_vir

class UCIS(CIS):
  def __init__(self, mol, verbose, save_states):
      CIS.__init__(self, mol, verbose, save_states)
      self.same_ab = False

  def _set_scf_reference(self):
      psi4.core.set_global_option('reference', 'uhf')

  def _set_beta(self, H, eri):
      self.Cb_occ = self.ref_wfn.Cb_subset("AO","OCC")
      self.Cb_vir = self.ref_wfn.Cb_subset("AO","VIR")
      self.Db = self.ref_wfn.Db()
      self.nbocc = self.ref_wfn.nbeta()
      self.nbvir = self.nmo - self.nbocc
      self.ndet = 1 + self.naocc * self.navir + self.nbocc * self.nbvir
      self.hamiltonian = numpy.zeros((self.ndet, self.ndet),numpy.float64)
      #
      Oa = self.Ca_occ.to_array(dense=True)
      Va = self.Ca_vir.to_array(dense=True)
      Ob = self.Cb_occ.to_array(dense=True)
      Vb = self.Cb_vir.to_array(dense=True)
      #
      #self.eri_ovov = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Ob, Vb, Ob, Vb, eri)
      #self.eri_OVov = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Va, Ob, Vb, eri)
      #self.eri_oovv = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Ob, Ob, Vb, Vb, eri)
      self.eri_ovov = four_index_transform(eri, Ob, Vb, Ob, Vb)
      self.eri_OVov = four_index_transform(eri, Oa, Va, Ob, Vb)
      self.eri_oovv = four_index_transform(eri, Ob, Ob, Vb, Vb)
      #
      #self.jk.C_clear()
      #self.jk.C_left_add(psi4.core.Matrix.from_array(self.Db, ""))
      #I = numpy.identity(self.Da.shape[0], numpy.float64)
      #self.jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      #self.jk.compute()
      #Jb= self.jk.J()[0].to_array(dense=True)
      #Kb= self.jk.K()[0].to_array(dense=True)
      ##
      #Gb = H + Ja + Jb - Kb
      Gb = self.ref_wfn.Fb().to_array(dense=True)
      #self.Fb_occ = numpy.einsum("ai,ab,bj->ij", Ob, Gb, Ob)
      #self.Fb_vir = numpy.einsum("ai,ab,bj->ij", Vb, Gb, Vb)
      self.Fb_occ = two_index_transform(Gb, Ob, Ob)
      self.Fb_vir = two_index_transform(Gb, Vb, Vb)


# - the same as above but in one class - #
class CIS_one_class:
  "Direct CIS Method"
  def __init__(self, mol, verbose=True, save_states=4):
      self.mol = mol
      self.save_states = save_states
      self.Ca_occ = None
      self.Cb_occ = None
      self.Ca_vir = None
      self.Cb_vir = None
      self.Da = None
      self.Db = None
      self.scf_e = None
      self.e_0 = None
      self.N = mol.nuclear_repulsion_energy()
      self.E = None
      self.W = None
      self.hamiltonian = None
      self.nmo = None
      self.naocc = None
      self.nbocc = None
      self.navir = None
      self.nbvir = None
      self.ndet = None
      self.jk = None
  def run(self):
      self._run_scf()
      self._prepare_for_cis()
      self._build_hamiltonian()
      self._diagonalize()
  def _prepare_for_cis(self):
      self.Ca_occ = self.ref_wfn.Ca_subset("AO","OCC")
      self.Cb_occ = self.ref_wfn.Cb_subset("AO","OCC")
      self.Ca_vir = self.ref_wfn.Ca_subset("AO","VIR")
      self.Cb_vir = self.ref_wfn.Cb_subset("AO","VIR")
      self.Da = self.ref_wfn.Da()
      self.Db = self.ref_wfn.Db()
      self.nmo = self.ref_wfn.nmo()
      self.naocc = self.ref_wfn.nalpha()
      self.nbocc = self.ref_wfn.nbeta()
      self.navir = self.nmo - self.naocc
      self.nbvir = self.nmo - self.nbocc
      self.ndet = 1 + self.naocc * self.navir + self.nbocc * self.nbvir
      self.hamiltonian = numpy.zeros((self.ndet, self.ndet),numpy.float64)
      #
      self.jk = psi4.core.JK.build(self.ref_wfn.basisset(), jk_type='direct')
      self.jk.set_memory(int(5e8))
      self.jk.initialize()
      #
      mints = psi4.core.MintsHelper(self.ref_wfn.basisset())
      eri = numpy.asarray(mints.ao_eri())
      Oa = self.Ca_occ.to_array(dense=True)
      Ob = self.Cb_occ.to_array(dense=True)
      Va = self.Ca_vir.to_array(dense=True)
      Vb = self.Cb_vir.to_array(dense=True)
      #
      self.eri_OVOV = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Va, Oa, Va, eri)
      self.eri_ovov = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Ob, Vb, Ob, Vb, eri)
      self.eri_OVov = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Va, Ob, Vb, eri)
      self.eri_OOVV = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Oa, Oa, Va, Va, eri)
      self.eri_oovv = numpy.einsum("ai,bj,ck,dl,abcd->ijkl", Ob, Ob, Vb, Vb, eri)
      #
      #H = self.ref_wfn.H().to_array(dense=True)
      #
      #self.jk.C_clear()
      #self.jk.C_left_add(psi4.core.Matrix.from_array(self.Da, ""))
      #I = numpy.identity(self.Da.shape[0], numpy.float64)
      #self.jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      #self.jk.compute()
      #Ja= self.jk.J()[0].to_array(dense=True)
      #Ka= self.jk.K()[0].to_array(dense=True)
      ##
      #self.jk.C_clear()
      #self.jk.C_left_add(psi4.core.Matrix.from_array(self.Db, ""))
      #self.jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      #self.jk.compute()
      #Jb= self.jk.J()[0].to_array(dense=True)
      #Kb= self.jk.K()[0].to_array(dense=True)
      #
      #Fa = H+Ja+Jb-Ka
      #Fb = Fa+Ka-Kb
      #print(Fa, self.ref_wfn.Fa().to_array(dense=True))
      Fa = self.ref_wfn.Fa().to_array(dense=True)
      Fb = self.ref_wfn.Fb().to_array(dense=True)
      self.Fa_occ = numpy.einsum("ai,ab,bj->ij", Oa, Fa, Oa)
      self.Fa_vir = numpy.einsum("ai,ab,bj->ij", Va, Fa, Va)
      self.Fb_occ = numpy.einsum("ai,ab,bj->ij", Ob, Fb, Ob)
      self.Fb_vir = numpy.einsum("ai,ab,bj->ij", Vb, Fb, Vb)
  def _run_scf(self):
      scf_e, wfn = psi4.energy('HF', molecule=self.mol, return_wfn=True)
      self.scf_e = scf_e
      self.e_0   = scf_e - self.N
      self.ref_wfn = wfn
  def _build_hamiltonian(self):
      # OO block
      self.hamiltonian[0, 0] = self.e_0
      # OS and SO blocks are zero
      # SS block
      off_a = self.naocc*self.navir
      off_b = self.nbocc*self.nbvir
      for i in range(self.naocc):
          for a in range(self.navir):
              ia = (self.navir)*i + a
              for j in range(self.naocc):
                  for b in range(self.navir):
                      jb = (self.navir)*j + b
                      v = 0.0
                      if (i==j) and (a==b): v+= self.e_0
                      if (i==j): v+= self.Fa_vir[a,b]
                      if (a==b): v-= self.Fa_occ[i,j]
                      v += self.eri_OVOV[i,a,j,b] - self.eri_OOVV[i,j,a,b] 
                      #
                      self.hamiltonian[1+ia,1+jb] = v
              for j in range(self.nbocc):
                  for b in range(self.nbvir):
                      jb = (self.nbvir)*j + b
                      v = self.eri_OVov[i,a,j,b]
                      self.hamiltonian[1+ia,1+jb+off_a] = v
                      self.hamiltonian[1+jb+off_a,1+ia] = v
      for i in range(self.nbocc):
          for a in range(self.nbvir):
              ia = (self.nbvir)*i + a
              for j in range(self.nbocc):
                  for b in range(self.nbvir):
                      jb = (self.nbvir)*j + b
                      v = 0.0
                      if (i==j) and (a==b): v+= self.e_0
                      if (i==j): v+= self.Fb_vir[a,b]
                      if (a==b): v-= self.Fb_occ[i,j]
                      v += self.eri_ovov[i,a,j,b] - self.eri_oovv[i,j,a,b] 
                      #
                      self.hamiltonian[1+ia+off_a,1+jb+off_a] = v
      del self.eri_ovov, self.eri_OVOV, self.eri_OVov
  def _diagonalize(self):
      E, W = numpy.linalg.eigh(self.hamiltonian)
      self.E = E
      self.W = W
