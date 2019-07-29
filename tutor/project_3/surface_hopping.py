#!/usr/bin/python3
"""Surface Hopping Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import math
import psi4
import numpy
from .trajectory import Aggregate, TimePoint, Trajectory
from .hamiltonian import Isolated_Hamiltonian

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
