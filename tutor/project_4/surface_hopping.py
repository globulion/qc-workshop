#!/usr/bin/python3
"""Surface Hopping Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from .aggregate import Aggregate
from .trajectory import TimePoint, Trajectory
from .hamiltonian import Isolated_Hamiltonian
from ..psithon.util import rearrange_eigenpairs, _reorder, check_sim

class System:
  def __init__(self, psi4_molecule, temperature=0.0, nstates=1, current_state=0):
      # Molecular aggregate
      self.aggregate = Aggregate(psi4_molecule)
      # Physical temperature
      self.temperature = temperature
      # Number of all electronic states
      self.nstates = nstates
      # Current electronic state index
      self.current_state = current_state
      # Hamiltonian of entire system
      self.hamiltonian = None
      
  def update(self, xyz):
      "Update geometry of all molecules in the system"
      self.aggregate.update(xyz)
      self.hamiltonian.computers[0].update(self.aggregate.qm.geometry())
      for i, computer in enumerate(self.hamiltonian.computers[1:]):
          computer.update(self.aggregate.bath[i].geometry())

  def set_hamiltonian(self, method_high, method_low=None):
      "Set the Hamiltonian of the system"
      if method_low is not None: raise NotImplementedError("Now only isolated Hamiltonian is implemented")
      self.hamiltonian = Isolated_Hamiltonian(method_high, self.aggregate, self.nstates)

class Units:
  fs2au = 4.1341373336493e+16 * 1.0e-15
  au2fs = 1./fs2au

class DynamicalSystem(System, Units):
  def __init__(self, aggregate, nstates=1, init_state=0, temperature=0.0, 
                     qm_method='myCIS', gs_method=None, dt_class=0.5, dt_quant=None, seed=0):
      System.__init__(self, aggregate, temperature, nstates, init_state)
      Units.__init__(self)
      numpy.random.seed(seed)

      # Trajectory of the system
      self.trajectory = Trajectory(self.aggregate)
      # Initial electronic state index
      self.init_state = init_state
      # Current electronic state index
      self.current_state = init_state
      # Dimension of quantum adiabatic space
      self.dim = nstates
      # Time steps
      self.dt_class = dt_class * self.fs2au
      self.dt_quant = dt_class * self.fs2au / 1000.0 if dt_quant is None else dt_quant * self.fs2au
      # Instantaneous total energy of the system
      self.energy = None
      # Instantaneous quantum amplitudes
      self.c = None
      # Instantaneous density matrix
      self.d = None

      # Set up the Hamiltonian      
      self.set_hamiltonian(qm_method)

  def run(self, nt, out_xvf='traj.dat', out_xyz='traj.xyz', center_mode='qm'):
      "Run the TSH dynamics"
      outf = open(out_xvf, 'w')
      outx = open(out_xyz, 'w')

      print(" Initial Conditions")
      self.set_initial_conditions()
      self.trajectory.save(outf)
      self.aggregate.save_xyz(outx, center_mode)

      for i in range(nt): 
          t = i*self.dt_class*self.au2fs
          print(" t = %13.3f [fs]  s = %2d  T = %6.1f [K]" % (t, self.current_state, self.temperature))
          self.propagate()
          self.trajectory.save(outf)
          self.aggregate.save_xyz(outx, center_mode)
          print(" Occupancies: ")
          print("%8.3f"*len(self.p) % tuple(self.d.real.diagonal()))
          print(" Transition Probabilities: ")
          print("%8.3f"*len(self.p) % tuple(self.p))

      outf.close()
      outx.close()

  def try_to_hop(self, g, e_curr, e_old): 
      "Try to hop from k to m"
      zeta = numpy.random.random()
      print("Try to hop! z = %7.4f" % zeta)
      e_k = e_old[self.current_state]
      m_hop = None
      sum_with_m = g.sum()
      for m in range(g.size):
          sum_without_m = sum_with_m - g[m]
          if sum_without_m < zeta <= sum_with_m:
             print("Almost hop!")
             e_m = e_curr[m]
             if e_k >= e_m:
                print("Hop %d --> %d" % (self.current_state, m))
                m_hop = m
             else: 
                print("Hop is frustrated!")
                m_hop = m
      if m_hop is not None:
         self.current_state = m_hop
         self._decoherence_correction()

  def set_initial_conditions(self):
      # amplitudes
      self.c = numpy.zeros(self.dim, dtype=numpy.complex128)
      self.c[self.init_state] = 1.0+0.0j
      # density matrix
      self._update_density_matrix(self.c)
      # Hamiltonian
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
      #self.trajectory.rescale_velocities(0.01)
      #
      self.energy = s.ci_e[self.init_state] + self.aggregate.all.nuclear_repulsion_energy() + self.trajectory.kinetic_energy()
      print("Temp init= ", self.trajectory.temperature())

  def propagate(self):
      dt = self.dt_class
      # [0] Grab current time point
      x_old = self.trajectory.point_last.x
      v_old = self.trajectory.point_last.v
      a_old = self.trajectory.point_last.f * self.aggregate._mrec[:, numpy.newaxis]
      s_old = self.trajectory.point_last.s
      c_old = self.c.copy()
      p_old = self.d.copy()
      e_el_k   = s_old.ci_e[self.current_state]
      Tkin_old = self.trajectory.kinetic_energy()
      nrep_old = self.aggregate.all.nuclear_repulsion_energy()

      # [6] Compute next positions
      x_new = x_old + dt*v_old + 0.5 * dt*dt * a_old
      self.update(psi4.core.Matrix.from_array(x_new))
      nrep_new = self.aggregate.all.nuclear_repulsion_energy()

      # [1] Compute next wavefunction
      self.hamiltonian.compute()
      s_new = self.hamiltonian.computers[0].ciwfn
      #f = self.hamiltonian.computers[0].forces
      #if 0:
      #   ci_e = s_new.ci_e                                                                                      
      #   ci_c_= s_new.ci_c
      #   ci_e, ci_c, sim = rearrange_eigenpairs(ci_c_, s_old.ci_c, ci_e, return_sim=True)
      #   s_new.ci_e = ci_e
      #   s_new.ci_c = ci_c
      #   log = check_sim(sim)
      #   if 'ERROR' in log: raise ValueError("Error in state rearrangement! %s\n\n%s" % (str(sim), str(ci_c_)))
      #   f = _reorder(f_new, sim)
      #   print(sim)

      # [2] Grab current Hamiltonian matrix in adiabatic basis
      H = numpy.diag(s_new.ci_e)
                                                                   
      # [3] Compute non-adiabatic coupling constants sigma_ij
      S = s_old.overlap(s_new)
      s =(S - S.T) / (2.0 * dt)
                                                                   
      # [4] Compute G operator matrix
      G = H - 1.0j*s.T
                                                                   
      # [5] Compute next dynamical amplitudes
      c_new = self._step_quantum(c_old, G)
      self._update_density_matrix(c_new)

      # [8] Compute probabilities from current state
      d_ii= self.d.diagonal()[self.current_state].real
      d_ij= self.d[self.current_state].real
      s_i = s[self.current_state]
      P = 2.0 * dt * d_ij * s_i / d_ii
      P[P<0.0] = 0.0
      norm = numpy.linalg.norm(P)
      #if norm>1.0: P/= norm
      self.p = P

      # [9] Hop if required
      self.try_to_hop(P, s_new.ci_e, s_old.ci_e)
      e_el_m = s_new.ci_e[self.current_state]

      # [7] Compute next forces and velocities
      f_new = self.hamiltonian.computers[0].forces[self.current_state]
      v_new = v_old + 0.5 * dt * (a_old + f_new * self.aggregate._mrec[:, numpy.newaxis])
      #v_new.fill(0.0)

      # [8] Save
      point_next = TimePoint(x_new, v_new, f_new, s_new)
      self.trajectory.add_point(point_next)
      e_kin_target = Tkin_old + e_el_k - e_el_m + nrep_old - nrep_new
      #print("T=",e_kin_target, e_el_k, e_el_m)
      if e_kin_target> 0.0: self.trajectory.rescale_velocities(e_kin_target)
      self.temperature = self.trajectory.temperature()


  def _density_matrix(self, c):
      return numpy.outer(self.c, self.c.conjugate())

  def _update_density_matrix(self, c):
      self.c = c.copy()
      self.d = self._density_matrix(self.c)

  def _step_quantum(self, c, G): 
      g, u = numpy.linalg.eig(-1.0j*G)
      e = numpy.linalg.multi_dot([u, numpy.diag(numpy.exp(g*self.dt_class)), u.T])
      return numpy.dot(c, e.T)

  def _decoherence_correction(self):
      self.c.fill(0.0+0.0j)
      self.c[self.current_state] = 1.0+0.0j
      self._update_density_matrix(self.c)

