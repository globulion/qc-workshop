#!/usr/bin/python3
"""
 Module for surface hopping method
"""
from abc import ABC, abstractmethod
import psi4
import numpy

class MCWavefunction:
  def __init__(self, ci_e, ci_c, ci_l, ca, cb, bfs): 
      self.ci_e = ci_e
      self.ci_c = ci_c
      self.ci_l = ci_l
      self.ca   = ca
      self.cb   = cb
      self.bfs  = bfs
      self.ndet = len(ci_l)
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
      S_nm  = numpy.zeros((self.ndet, other.det), numpy.float64)
      # sum over all Slater determinants
      for n in range(self.ndet):
          cai, cbi = self._construct_c(n)
          for m in range(other.ndet):
              caj, cbj = other._construct_c(m)
              D_ij_a = numpy.linalg.multi_dot(cai.T, sab, caj)
              D_ij_b = numpy.linalg.multi_dot(cbi.T, sab, cbj)
              S[n,m] = numpy.linalg.det(D_ij_a) * numpy.linalg.det(D_ij_b)
      return S_nm
  def _construct_c(self, n):
      o, v = self.ci_l[n]
      ca   = self.ca.copy()
      cb   = self.cb.copy()
          
      return ca, cb

class TimePoint:
  def __init__(self, x, v, f, s=None):
      self.x = x
      self.v = v
      self.f = f
      self.s = s
  def save(self, out, n=None):
      if n is not None:
         x = self.x[n]
         v = self.v[n]
         f = self.f[n]
      else:
         x = self.x
         v = self.v
         f = self.f
      log = "%13.6f"*3 % tuple(x) + '\n'
      log+= "%13.6f"*3 % tuple(v) + '\n'
      log+= "%13.6f"*3 % tuple(f) + '\n'
      out.write(log)

class Trajectory(ABC):
  def __init__(self, natoms, nstates):
      ABC.__init__(self)
      self.natoms = natoms
      self.nstates= nstates
      self.point_last = None
      self.point_prev = None
  def add_point(self, point):
      self.point_prev = self.point_last
      self.point_last = point
  def save(self, out, n=None): self.point_last.save(out, n)
      
  @abstractmethod
  def x(self, i, n=None): pass
  @abstractmethod
  def v(self, i, n=None): pass
  @abstractmethod
  def f(self, i, n=None): pass

class QuantumTrajectory(Trajectory):
  def __init__(self, natoms, nstates):
      Trajectory.__init__(self, natoms, nstates)
  def x(self, i, n): return self.trajectory[i].x[n]
  def v(self, i, n): return self.trajectory[i].v[n]
  def f(self, i, n): return self.trajectory[i].f[n]

def ClassicalTrajectory(Trajectory):
  def __init__(self, natoms): 
      Trajectory.__init__(self, natoms, 1)
  def x(self, i): return self.trajectory[i].x
  def v(self, i): return self.trajectory[i].v
  def f(self, i): return self.trajectory[i].f


class System:
  def __init__(self, aggregate, temp=0.0):
      self.aggregate = aggregate
      self.temp = temp
      self.hamiltonian = None
      self.molecule_qm = None
      self.molecule_bath = []
      self.trajectories_qm = []
      self.trajectories_bath = []
      self.update_aggregate(aggregate.geometry())

  def hop(self, nstate): pass

  def canonicalize_velocities(self, temp): pass

  def update_aggregate(self, xyz):
      self.aggregate.set_geometry(xyz)
      self.molecule_qm = self.aggregate.extract_subset(1)
      self.molecule_bath = []
      for i in range(1, self.n_molecules):
          self.molecule_bath.append(self.n_molecules.extract_subset(i))

  def set_hamiltonian(self, qm_method):
      self.hamiltonian = Isolated_Hamiltonian(qm_method, self.molecule_qm)


class Computer(ABC):
  def __init__(self, molecule):
      ABC.__init__(self)
      #
      self.molecule = molecule
      self.nuclear_repulsion_energy = molecule.nuclear_repulsion_energy()
      #
      self.state_amplitudes = None
      self.state_energies = None
      self.state_forces = None

  @classmethod
  def create(cls, method, molecule):
      m = method.lower()
      if m == "psi4scf": return psi4SCF_Computer(molecule)
      elif m == "mycis": return myCIS_Computer(molecule)
      else: raise ValueError("Wrong method chosen for computer")

  def update(self, xyz):#OK
      self.molecule.set_geometry(xyz)
      self.molecule.update_geometry()

  def compute(self): 
      self.state_electronic_energies, self.state_amplitudes = self._compute_energy()
      self.state_forces = self._compute_forces()

  @abstractmethod
  def _compute_energy(self): pass

  def _compute_forces(self):#OK
      xyz = self.molecule.geometry()
      d = 0.000001
      f = []
      e_0 = self.state_electronic_energies + self.nuclear_repulsion_energy
      n_states = len(e_0)
      for a in range(self.molecule.natom()):
       for x in range(3):
          xyz_xa1 = xyz.clone()
          xyz_xa1.set(a, x, xyz.get(a, x) + d)
          #
          self.update(xyz_xa1)
          e_xa1, w_xa1 = self._compute_energy()
          f_xa1 = (e_xa1 + self.molecule.nuclear_repulsion_energy() - e_0) / d
          #
          f.append(f_xa1)
      self.update(xyz)
      f = numpy.array(f).reshape(self.molecule.natom(),3,n_states).transpose(2,0,1)
      return f



class psiSCF_Computer(Computer):
  def __init__(self, molecule): 
      Computer.__init__(self, molecule)
  def _compute_energy(self):
      W = numpy.array([1.0])
      e = psi4.energy("SCF", molecule=self.molecule, return_wfn=False)
      E = numpy.array([e - self.molecule.nuclear_repulsion_energy()])
      return E, W

class CIS_Computer(Computer):
  def __init__(self, molecule):
      Computer.__init__(self, molecule) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule):
      CIS_Computer.__init__(self, molecule)
  def _compute_energy(self):
      import lib.cis
      cis = lib.cis.CIS(self.molecule, verbose=False, save_states=4)
      return cis.E, cis.W

class Hamiltonian(ABC):
  def __init__(self, method, molecule):
      ABC.__init__(self)
      self.computer = Computer.create(method, molecule)
      self.H = None
      self.F = None

  @classmethod
  def create(cls, molecule, method): pass

  def compute(self, bath_hamiltonians=None):
      "Solve Schrodinger Equation"
      # unperturbed Hamiltonian
      self.computer.compute()
      # environment
      if bath_hamiltonians is not None:
         self.iterate(h, bath_hamiltonians) 
      return
 
  @abstractmethod
  def iterate(self, h_0, h): pass

class Isolated_Hamiltonian(Hamiltonian):
  def __init__(self, molecule, method):
      Hamiltonian.__init__(self, molecule, method)
  def iterate(self, h_0, h): raise NotImplementedError("This Hamiltonian is isolated")

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

class Dynamics:
  def __init__(self, aggregate, qm_method='myCIS', gs_method=None, dt_class=0.00005, dt_quant=None,
                     init_state=0, max_states=4):
      self.dim = max_states
      self.dt_class = dt_class
      if dt_quant is None: dt_quant = dt_class
      self.system = System(aggregate)
      self.system.set_hamiltonian(qm_method)
      self.ct = numpy.zeros(max_states, dtype=numpy.complex128)
      self.ct[init_state] = 1.0

  def run(self, time):
      self._run(time)

  def _run(self, time): pass
      
  def _step_classical(self): pass

  def _step_quantum(self, c, G): 
      g, u = numpy.linalg.eig(G)
      e = numpy.linalg.multi_dot([u, numpy.exp(-1.0j*g*self.dt_class), u.T])
      return numpy.dot(c, e)
