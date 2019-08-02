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
  def save(self, out):#TODO
      raise NotImplementedError

class Trajectory:
  def __init__(self, aggregate):
      self.aggregate = aggregate
      self.natoms = aggregate.all.natom()
      self.point_last = None
      self.point_prev = None

  def add_point(self, point):#TODO
      raise NotImplementedError

  def set_points(self, prev, last):#TODO
      raise NotImplementedError

  def save(self, out): self.point_last.save(out)

  def kinetic_energy(self):#TODO
      "Compute kinetic energy of nuclei"
      raise NotImplementedError

  def rescale_velocities(self, e_kin):#TODO
      "Berendsen v-rescale thermostat"
      raise NotImplementedError

  def canonicalize_velocities(self, temp):#TODO
      raise NotImplementedError

  def temperature(self):#TODO
      raise NotImplementedError

