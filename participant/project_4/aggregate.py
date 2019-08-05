#!/usr/bin/python3
"""Aggregate Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
import psi4, numpy

class Aggregate:
  def __init__(self, psi4_molecule):
      self.all = psi4_molecule
      self.qm = psi4_molecule.extract_subsets(1)
      self.nfrags = psi4_molecule.nfragments()
      self.bath = [] if self.nfrags == 1 else [psi4_molecule.extract_subsets(2+i) for i in range(self.nfrags-1)]
      self._mrec = 1./(numpy.array([self.all.mass(i) for i in range(self.all.natom())])) * psi4.constants.au2amu
  def update(self, xyz):#TODO
      raise NotImplementedError
  def save_xyz(self, out, center_mode='qm'):#TODO
      raise NotImplementedError

