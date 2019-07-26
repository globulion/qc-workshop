#!/usr/bin/python3
"""
 Module for surface hopping method
"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from ..project_2.cis import MCState, CIS

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
      self.state_electronic_energies = None
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
      self.nuclear_repulsion_energy = self.molecule.nuclear_repulsion_energy()

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
      cis = CIS.create(self.molecule, verbose=False, save_states=4, reference='rhf')
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


class TimePoint:
  def __init__(self, x, v, f, s=None):
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
  def rescale_velocities(self, e_kin):
      "Berendsen v-rescale thermostat"
      e = 0.0
      for i in range(self.natoms):
          vi = self.point_last.v[i]
          e += self.molecules.mass(i) * numpy.dot(vi, vi) / psi4.constants.au2amu
      e/= 2.0 
      alpha = math.sqrt(e_kin/e)
      self.point_last.v*= alpha

class Aggregate:
  def __init__(self, molecule):
      self.all = molecule
      self.qm = molecule.extract_subsets(1)
      self.nfrags = molecule.nfragments()
      self.bath = [] if self.nfrags == 1 else [molecule.extract_subsets(2+i) for i in range(self.nfrags-1)]
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
      
  def canonicalize_velocities(self, temp): 
      print("Warning: no Maxwell-Boltzmann distribution is applied yet")

  def rescale_velocities(self, e_kin): 
      self.trajectory.rescale_velocities(e_kin)

  def update_aggregate(self, xyz):
      self.aggregate.update(xyz)
      self.hamiltonian.computer.update(self.aggregate.qm.geometry())

  def set_hamiltonian(self, qm_method):
      self.hamiltonian = Isolated_Hamiltonian(qm_method, self.aggregate.qm)

class Units:
  fs2au = 4.1341373336493e+16 * 1.0e-15
  au2fs = 1./fs2au

class DynamicalSystem(System, Units):
  def __init__(self, aggregate, temperature=0.0, nstates=1,
                     qm_method='myCIS', gs_method=None, dt_class=0.5, dt_quant=None,
                     init_state=0, max_states=5):
      System.__init__(self, aggregate, temperature, nstates)
      Units.__init__(self)
      #
      self.dim = max_states
      self.init_state = init_state
      self.current_state = init_state
      self.temperature = temperature
      self.dt_class = dt_class * self.fs2au
      if dt_quant is None: dt_quant = dt_class
      self.set_hamiltonian(qm_method)
      #
      self.c = None
      self.d = None

  def run(self, nt, time=None, out='traj.dat', out_xyz='traj.xyz'):
      outf = open(out, 'w')
      outx = open(out_xyz, 'w')
      print(" Initial Conditions")
      self.set_initial_conditions()
      self.trajectory.save(outf)
      self.aggregate.save_xyz(outx)

      #nt = time/self.dt_class
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
          if g[:(state+1)].sum() < zeta < g[:(state+2)].sum():
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
      v = numpy.random.random((self.aggregate.all.natom(), 3))
      f = self.hamiltonian.computer.state_forces[self.init_state]
      s = self.hamiltonian.computer.states
      point_init = TimePoint(x, v, f, s)
      #
      self.trajectory.add_point(point_init)
      self.canonicalize_velocities(self.temperature)
      self.rescale_velocities(0.0)

  def propagate(self):
      dt = self.dt_class
      # [0] Grab current time point
      x_old = self.trajectory.point_last.x
      v_old = self.trajectory.point_last.v
      a_old = self.trajectory.point_last.f * self._m[:, numpy.newaxis]
      s_old = self.trajectory.point_last.s
      w_old = MCState(s_old)
      c_old = self.c.copy()
      p_old = self.d.copy()

      # [1] Compute next x
      x_new = x_old + dt*v_old + 0.5 * dt*dt * a_old
      self.update_aggregate(psi4.core.Matrix.from_array(x_new))

      # [2] Compute next wavefunction and forces
      self.hamiltonian.compute()
      s_new = self.hamiltonian.computer.states
      f_new = self.hamiltonian.computer.state_forces[self.current_state]

      # [3] Compute next velocities
      v_new = v_old + 0.5 * dt * (a_old + f_new * self._m[:, numpy.newaxis])
      v_new.fill(0.0)

      # [4] Grab previous Hamiltonian matrix in adiabatic basis
      H = numpy.diag(self.hamiltonian.computer.state_electronic_energies)

      # [4] Compute non-adiabatic coupling constants
      w_new = MCState(s_new) 
      S = w_old.overlap(w_new)
      s =(S - S.T) / (2.0 * dt)

      # [5] Compute G operator matrix
      G = H - 1.0j*s

      # [6] Compute next dynamical amplitudes
      c_new = self._step_quantum(c_old, G)
      self.c = c_new.copy()
      self._update_density_matrix()
      p_new = self.d

      # [7] Compute probabilities
      P = 2.0 * dt * p_new.real * s / p_new.diagonal()
      P = P[self.current_state]
      P[P<0.0] = 0.0

      # [8] Hop if required
      self.hop(P)

      # [9] Save
      point_next = TimePoint(x_new, v_new, f_new, s_new)
      self.trajectory.add_point(point_next)


  def _density_matrix(self, c):
      return numpy.outer(self.c, self.c.conjugate())

  def _update_density_matrix(self):
      self.d = self._density_matrix(self.c)

  def _step_quantum(self, c, G): 
      g, u = numpy.linalg.eig(G)
      e = numpy.linalg.multi_dot([u, numpy.diag(numpy.exp(-1.0j*g*self.dt_class)), u.T])
      return numpy.dot(c, e)
