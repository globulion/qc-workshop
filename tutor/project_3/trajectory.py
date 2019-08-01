#!/usr/bin/python3
"""Trajectory Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
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
  def __init__(self, aggregate):
      self.aggregate = aggregate
      self.natoms = aggregate.all.natom()
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
          e += self.aggregate.all.mass(i) * numpy.dot(vi, vi) / psi4.constants.au2amu
      e/= 2.0 
      return e
  def rescale_velocities(self, e_kin):
      "Berendsen v-rescale thermostat"
      assert(e_kin >= 0.0)
      e = self.kinetic_energy()
      #print(e_kin, e)
      if e> 0.0: alpha = math.sqrt(e_kin/e)
      else: alpha = 0.0
      #print("Alpha=",alpha)
      self.point_last.v*= alpha
  def canonicalize_velocities(self, temp): 
      t = self.temperature()
      alpha = math.sqrt(temp/t)
      self.point_last.v*= alpha
      print("Warning: no Maxwell-Boltzmann distribution is applied yet")
  def temperature(self):
      kb = psi4.constants.kb / psi4.constants.hartree2J
      t = 2.0/(3.0*self.natoms) * self.kinetic_energy() / kb
      #print("conversion", 2.0/(3.0*self.natoms)/kb)
      return t
