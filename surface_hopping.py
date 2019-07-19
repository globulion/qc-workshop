#!/usr/bin/python3
"""
 Module for surface hopping method
"""
class Trajectory:
  def __init__(self):

class System:
  def __init__(self, aggregate, temp=0.0):
      self.aggregate = aggregate
      self.temp = temp

  def hop(self, nstate): 

  def canonicalize_velocities(self, temp):

class Hamiltonian(ABC):
  def __init__(self, aggregate):

  @staticmethod
  def create():

  def update_aggregate(self, xyz):


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
  def __init__(self, aggregate, dt_class=0.00005, dt_quant=None):
      self.dt_class = dt_class
      if dt_quant is None: dt_quant = dt_class
      self.system = System(aggregate)

  def _run(self, time):
      
  def _step_classical(self):

  def _step_quantum(self):
