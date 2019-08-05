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
  def update(self, xyz):
      self.all.set_geometry(xyz)
      self.qm = self.all.extract_subsets(1)
      self.bath = [self.all.extract_subsets(2+i) for i in range(self.nfrags-1)]
  def save_xyz(self, out, center_mode='qm'):
      geom = self.all.geometry()
      geom.scale(psi4.constants.bohr2angstroms)
      if   center_mode is None         : com = [0.0,0.0,0.0]
      elif center_mode.lower() == 'qm' : com = self.qm.center_of_mass()
      elif center_mode.lower() == 'all': com = self.all.center_of_mass()
      elif isinstance(center_mode, int): com = self.bath[center_mode].center_of_mass()
      else: raise ValueError("Centering mode - %s - is not supported" % center_mode)
      out.write("%d\n\n" % self.all.natom())
      for i in range(self.all.natom()):                                                              
          sym = self.all.label(i)
          out.write("%s %10.6f %10.6f %10.6f\n"%(sym,geom.get(i,0)-com[0],geom.get(i,1)-com[1],geom.get(i,2)-com[2]))
