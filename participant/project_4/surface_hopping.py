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

  def try_to_hop(self, g, e_curr, e_old):#TODO
      "Try to hop from k to m"
      raise NotImplementedError

  def set_initial_conditions(self):#TODO
      raise NotImplementedError

  def propagate(self):#TODO
      raise NotImplementedError

  def _density_matrix(self, c):#TODO
      raise NotImplementedError

  def _update_density_matrix(self, c):#TODO
      raise NotImplementedError

  def _step_quantum(self, c, G):#TODO
      raise NotImplementedError

  def _decoherence_correction(self):#TODO
      raise NotImplementedError

