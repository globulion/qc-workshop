#!/usr/bin/python3
"""Quantum Computer Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
from abc import ABC, abstractmethod
import psi4
import math
import numpy

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

class Trajectory:
  def __init__(self, molecules):
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
