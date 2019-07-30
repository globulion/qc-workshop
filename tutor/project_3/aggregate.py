#!/usr/bin/python3
"""Aggregate Module.

.. moduleauthor:: Bartosz Błasiak <blasiak.bartosz@gmail.com>

"""
import psi4

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
