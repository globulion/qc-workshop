#!/usr/bin/python3
"""Surface Hopping Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from ..project_2.cis import CIS, HF_CIWavefunction, CIS_CIWavefunction

__all__ = ["SlaterDeterminant",
              "Reference_SlaterDeterminant",
              "Single_SlaterDeterminant",
           "CIWavefunction",
              "HF_CIWavefunction",
              "CIS_CIWavefunction"]

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


class Computer(ABC):
  def __init__(self, molecule):
      ABC.__init__(self)
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
      self.molecule.update_geometry()
      self.nuclear_repulsion_energy = self.molecule.nuclear_repulsion_energy()

  def compute(self): 
      self.ciwfn = self._compute_energy()
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
      d = 0.000001
      f = []
      e_0 = self.ciwfn.ci_e + self.nuclear_repulsion_energy
      n_states = len(e_0)
      for a in range(self.molecule.natom()):
          for x in range(3):                                                           
             xyz_xa1 = xyz.clone()
             xyz_xa1.set(a, x, xyz.get(a, x) + d)
             #
             self.update(xyz_xa1)
             ciwfn = self._compute_energy()
             f_xa1 = (ciwfn.ci_e + self.molecule.nuclear_repulsion_energy() - e_0) / d
             #
             f.append(f_xa1)
      self.update(xyz)
      f = -numpy.array(f).reshape(self.molecule.natom(),3,n_states).transpose(2,0,1)
      return f


class CIS_Computer(ExcitedState_Computer):
  def __init__(self, molecule, nstates):
      ExcitedState_Computer.__init__(self, molecule, nstates) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule, nstates):
      CIS_Computer.__init__(self, molecule, nstates)
  def _compute_energy(self):
      cis = CIS.create(self.molecule, verbose=False, save_states=self.nstates-1, reference='uhf')
      cis.run()
      state = CIS_CIWavefunction(cis.ref_wfn, cis.E, cis.W)
      return state

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


class TimePoint:
  def __init__(self, x, v, f, s):
      self.x = x.copy()
      self.v = v.copy()
      self.f = f.copy()
      self.s = s
      self.time_id = 0
  def save(self, out):
      self.time_id+= 1
      log = "%3d" % self.time_id
      n = self.x.size
      #log+= "%13.6f"*n % tuple(self.x.ravel()) + '\n'
      #log+= "%13.6f"*n % tuple(self.v.ravel()) + '\n'
      log+= "%13.6f"*3 % tuple(numpy.linalg.norm(self.f, axis=1)) + '\n'
      out.write(log)

class Trajectory(ABC):
  def __init__(self, molecules):
      ABC.__init__(self)
      self.molecules = molecules
      self.natoms = molecules.natom()
      self.point_last = None
      self.point_prev = None
  def add_point(self, point):
      self.point_prev = self.point_last
      self.point_last = point
  def set_points(self, prev, last):
      self.point_prev = prev 
      self.point_last = last
  def save(self, out): self.point_last.save(out)
  def kinetic_energy(self):
      "Compute kinetic energy of nuclei"
      e = 0.0
      for i in range(self.natoms):
          vi = self.point_last.v[i]
          e += self.molecules.mass(i) * numpy.dot(vi, vi) / psi4.constants.au2amu
      e/= 2.0 
      return e
  def rescale_velocities(self, e_kin):
      "Berendsen v-rescale thermostat"
      #if e_kin < 0: e_kin = 0.0
      e = self.kinetic_energy()
      if e> 0.0: alpha = math.sqrt(e_kin/e)
      else: alpha = 0.0
      self.point_last.v*= alpha
  def canonicalize_velocities(self, temp): 
      e_kin = self.kinetic_energy()
      t = 2.0/3.0 * e_kin / psi4.constants.kb * psi4.constants.hartree2J
      alpha = math.sqrt(temp/t)
      self.point_last.v*= alpha
      print("Warning: no Maxwell-Boltzmann distribution is applied yet")

class Aggregate:
  def __init__(self, psi4_molecule):
      self.all = psi4_molecule
      self.qm = psi4_molecule.extract_subsets(1)
      self.nfrags = psi4_molecule.nfragments()
      self.bath = [] if self.nfrags == 1 else [psi4_molecule.extract_subsets(2+i) for i in range(self.nfrags-1)]
  def update(self, xyz):
      self.all.set_geometry(xyz)
      self.qm = self.all.extract_subsets(1)
      self.bath = [self.all.extract_subsets(2+i) for i in range(self.nfrags-1)]
  def save_xyz(self, out):
      geom = self.all.geometry()
      geom.scale(psi4.constants.bohr2angstroms)
      out.write("%d\n\n" % self.all.natom())
      for i in range(self.all.natom()):                                                              
          sym = self.all.label(i)
          out.write("%s %10.6f %10.6f %10.6f\n"%(sym,geom.get(i,0),geom.get(i,1),geom.get(i,2)))
     
class System:
  def __init__(self, molecule, temperature=0.0, nstates=1):
      self.aggregate = Aggregate(molecule)
      self.temperature = temperature
      self.nstates = nstates
      self.current_state = None
      self.hamiltonian = None
      #
      self.trajectory = Trajectory(self.aggregate.all)
      #
      self._m = 1./(numpy.array([self.aggregate.all.mass(i) for i in range(self.aggregate.all.natom())])) * psi4.constants.au2amu
      
  def update(self, xyz):
      self.aggregate.update(xyz)
      self.hamiltonian.computers[0].update(self.aggregate.qm.geometry())
      for i, computer in enumerate(self.hamiltonian.computers[1:]):
          computer.update(self.aggregate.bath[i].geometry())

  def set_hamiltonian(self, method_high, method_low=None):
      self.hamiltonian = Isolated_Hamiltonian(method_high, self.aggregate, self.nstates)

class Units:
  fs2au = 4.1341373336493e+16 * 1.0e-15
  au2fs = 1./fs2au

class DynamicalSystem(System, Units):
  def __init__(self, aggregate, nstates=1, init_state=0, temperature=0.0, 
                     qm_method='myCIS', gs_method=None, dt_class=0.5, dt_quant=None, seed=0):
      System.__init__(self, aggregate, temperature, nstates)
      Units.__init__(self)
      numpy.random.seed(seed)
      #
      self.init_state = init_state
      self.current_state = init_state
      self.temperature = temperature
      self.dt_class = dt_class * self.fs2au
      if dt_quant is None: dt_quant = dt_class
      self.set_hamiltonian(qm_method)
      self.dim = self.hamiltonian.computers[0].nstates
      #
      self.c = None
      self.d = None
      #
      self.energy = None

  def run(self, nt, out='traj.dat', out_xyz='traj.xyz'):
      outf = open(out, 'w')
      outx = open(out_xyz, 'w')

      print(" Initial Conditions")
      self.set_initial_conditions()
      self.trajectory.save(outf)
      self.aggregate.save_xyz(outx)

      for i in range(nt): 
          t = i*self.dt_class*self.au2fs
          print(" t = %13.3f [fs]  s = %2d" % (t, self.current_state))
          self.propagate()
          self.trajectory.save(outf)
          self.aggregate.save_xyz(outx)

      outf.close()
      outx.close()

  def hop(self, g): 
      zeta = numpy.random.random()
      for state in range(g.size):
          if g[:(state)].sum() < zeta < g[:(state+1)].sum():
             self.current_state = state

  def set_initial_conditions(self):
      # amplitudes
      self.c = numpy.zeros(self.dim, dtype=numpy.complex128)
      self.c[self.init_state] = 1.0+0.0j
      # density matrix
      self._update_density_matrix()
      # hamiltonian
      self.hamiltonian.compute()
      # phase space
      x = self.aggregate.all.geometry().to_array(dense=True)
      v = numpy.random.random((self.aggregate.all.natom(), 3)) - 1.0
      f = self.hamiltonian.computers[0].forces[self.init_state]
      s = self.hamiltonian.computers[0].ciwfn
      point_init = TimePoint(x, v, f, s)
      #
      self.trajectory.add_point(point_init)
      self.trajectory.canonicalize_velocities(self.temperature)
      self.trajectory.rescale_velocities(0.01)
      #
      self.energy = s.ci_e[self.init_state] + self.aggregate.all.nuclear_repulsion_energy() + self.trajectory.kinetic_energy()

  def propagate(self):
      dt = self.dt_class
      # [0] Grab current time point
      x_old = self.trajectory.point_last.x
      v_old = self.trajectory.point_last.v
      a_old = self.trajectory.point_last.f * self._m[:, numpy.newaxis]
      s_old = self.trajectory.point_last.s
      c_old = self.c.copy()
      p_old = self.d.copy()

      # [1] Compute next x
      x_new = x_old + dt*v_old + 0.5 * dt*dt * a_old
      self.update(psi4.core.Matrix.from_array(x_new))

      # [2] Compute next wavefunction and forces
      self.hamiltonian.compute()
      s_new = self.hamiltonian.computers[0].ciwfn
      f_new = self.hamiltonian.computers[0].forces[self.current_state]

      # [3] Compute next velocities
      v_new = v_old + 0.5 * dt * (a_old + f_new * self._m[:, numpy.newaxis])

      # [4] Grab previous Hamiltonian matrix in adiabatic basis
      H = numpy.diag(s_old.ci_e)

      # [4] Compute non-adiabatic coupling constants sigma_ij
      S = s_old.overlap(s_new)
      s =(S - S.T) / (2.0 * dt)

      # [5] Compute G operator matrix
      G = H - 1.0j*s.T

      # [6] Compute next dynamical amplitudes
      c_new = self._step_quantum(c_old, G)
      self.c = c_new.copy()
      self._update_density_matrix()
      p_new = self.d

      # [7] Compute probabilities from current state
      print(self.d.real.diagonal())
      d_ii= self.d.real.diagonal()[self.current_state]
      d_ij= self.d[self.current_state].real
      s_i = s[self.current_state]
      P = 2.0 * dt * d_ij * s_i / d_ii
      P[P<0.0] = 0.0
      print(P)

      # [8] Hop if required
      point_next = TimePoint(x_new, v_new, f_new, s_new)
      self.trajectory.add_point(point_next)
      e_kin_target = self.energy - self.aggregate.all.nuclear_repulsion_energy() - s_new.ci_e[self.current_state]
      if e_kin_target > 0: 
         self.hop(P)
         self.trajectory.rescale_velocities(e_kin_target)


  def _density_matrix(self, c):
      return numpy.outer(self.c, self.c.conjugate())

  def _update_density_matrix(self):
      self.d = self._density_matrix(self.c)

  def _step_quantum(self, c, G): 
      g, u = numpy.linalg.eig(G)
      e = numpy.linalg.multi_dot([u, numpy.diag(numpy.exp(-1.0j*g*self.dt_class)), u.T])
      return numpy.dot(c, e)
