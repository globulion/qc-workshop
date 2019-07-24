#!/usr/bin/python3
"""
 Module for surface hopping method
"""
from abc import ABC, abstractmethod
import psi4
import numpy


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

class PhaseSpace(ABC):
  def __init__(self, natoms, nstates):
      ABC.__init__(self)
      self.natoms = natoms
      self.nstates = nstates
  def set_x(self, x): self.x = x
  def set_v(self, v): self.v = v
  def set_f(self, f): self.f = f

  def rescale_velocities(self, temp): pass

class Quantum_PhaseSpace(PhaseSpace):
  def __init__(self, natoms, nstates):
      PhaseSpace.__init__(self, natoms, nstates) 
      self.x = numpy.zeros((self.nstates, self.natoms, 3), dtype=numpy.float64)
      self.v = numpy.zeros((self.nstates, self.natoms, 3), dtype=numpy.float64)
      self.f = numpy.zeros((self.nstates, self.natoms, 3), dtype=numpy.float64)
  def set_xn(self, x, n): self.x[n] = x.copy()
  def set_vn(self, v, n): self.v[n] = v.copy()
  def set_fn(self, f, n): self.f[n] = f.copy()

class Classical_PhaseSpace(PhaseSpace):
  def __init__(self, natoms):
      PhaseSpace.__init__(self, natoms, 1) 
      self.x = numpy.zeros((self.natoms, 3), dtype=numpy.float64)
      self.v = numpy.zeros((self.natoms, 3), dtype=numpy.float64)
      self.f = numpy.zeros((self.natoms, 3), dtype=numpy.float64)
     
class System:
  def __init__(self, aggregate, temp=0.0, nstates=1):
      self.aggregate = aggregate
      self.temp = temp
      self.nstates = 1
      self.hamiltonian = None
      #
      self.update_aggregate(aggregate.geometry())
      self.molecule_qm = aggregate.extract_subsets(1)
      self.molecule_bath = [] if aggregate.nfrag() == 1 else [aggregate.extract_subsets(2+i) for i in range(aggregate.nfrag()-1)]
      self.trajectories_qm = []
      self.trajectory_bath = []
      self.active_state = None
      self.phase_space_qm = Quantum_PhaseSpace(self.molecule_qm.natom(), nstates)
      self.phase_space_cl = None
      

  def hop(self, nstate): 
      xyz = self.aggregate.geometry().to_array(dense=True)
      xyz[:self.molecule_qm.natom(),:] = self.phase_space_qm.x[nstate]
      self.aggregate.update_aggregate(psi4.core.Matrix.from_array(xyz))

  def canonicalize_velocities(self, temp): 
      print("Warning: no Maxwell-Boltzmann distribution is applied yet")
  def rescale_velocities(self, temp): 
      self.phase_space_qm.rescale_velocities(temp)
      if self.phase_space_cl is not None: self.phase_space_cl.rescale_velocities(temp)

  def update_aggregate(self, xyz):
      self.aggregate.set_geometry(xyz)
      self.molecule_qm = self.aggregate.extract_subset(1)
      self.molecule_bath = []
      for i in range(1, self.n_molecules):
          self.molecule_bath.append(self.n_molecules.extract_subset(i))

  def set_hamiltonian(self, qm_method):
      self.hamiltonian = Isolated_Hamiltonian(qm_method, self.molecule_qm)

  def propagate(self, dt, init=False):
      if init:
         # quantum system
         x = self.aggregate.extract_subsets(1).geometry().to_array(dense=True)
      


class Computer(ABC):
  def __init__(self, molecule):
      ABC.__init__(self)
      #
      self.molecule = molecule
      self.nuclear_repulsion_energy = molecule.nuclear_repulsion_energy()
      #
      self.states = None
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

class ExcitedState_Computer(Computer):
  def __init__(self, molecule):
      Computer.__init__(self, molecule)

class CIS_Computer(ExcitedState_Computer):
  def __init__(self, molecule):
      ExcitedState_Computer.__init__(self, molecule) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule):
      CIS_Computer.__init__(self, molecule)
  def _compute_energy(self):
      import lib.cis
      cis = lib.cis.CIS.create(self.molecule, verbose=False, save_states=4, reference='rhf')
      cis.run()
      self.states = cis
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
      "Solve Schrodinger Equation for Extended Molecular Aggregate"
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
                     init_state=0, max_states=4, temperature=0.0):
      self.dim = max_states
      self.init_state = init_state
      self.temperature = temperature
      self.dt_class = dt_class
      if dt_quant is None: dt_quant = dt_class
      self.system = System(aggregate)
      self.system.set_hamiltonian(qm_method)
      self.state_current  = None
      self.state_previous = None

  def run(self, time):
      self._run(time)
  def _run(self, time): 
      self._set_initial_conditions()
   
      nt = time/self.dt_class
      for i in range(nt): pass

  def _set_initial_conditions(self):
      # amplitudes
      self.c = numpy.zeros(self.dim, dtype=numpy.complex128)
      self.c[self.init_state] = 1.0+0.0j
      # density matrix
      self._update_density_matrix()
      # velocities
      v = numpy.random.random((self.system.aggregate.natom(), 3))
      self.system.aggregate.set_geometry(psi4.core.Matrix.from_array(v))
      self.system.canonicalize_velocities(self.temperature)
      # hamiltonian
      self.system.hamiltonian.compute()
      self.state_previous = self.system.hamiltonian.computer.states
      #
      self.system.propagate(self.dt_class, init=True)
      self.state_current = self.system.hamiltonian.computer.states
      # 
      
  def _step_classical(self): pass

  def _compute_G(self): pass
  def _step_quantum(self, c, G): 
      g, u = numpy.linalg.eig(G)
      e = numpy.linalg.multi_dot([u, numpy.exp(-1.0j*g*self.dt_class), u.T])
      return numpy.dot(c, e)

  def _density_matrix(self, c):
      return numpy.outer(self.c, self.c.conjugate())
  def _update_density_matrix(self):
      self.d = self._density_matrix(self.c)

