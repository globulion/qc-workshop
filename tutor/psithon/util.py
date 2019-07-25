#!/usr/bin/python3
"""
 Workshop Utilities
"""
import numpy
import psi4

def matrix_power(P, a):
    "Power of a real symmetric matrix P^a"
    p, U = numpy.linalg.eigh(P)
    #p = abs(p)
    p[p<0.0] = 0.0
    #if (p<=0.0).any(): raise ValueError(" Matrix must be positive-definite!")
    Pa = numpy.linalg.multi_dot([U, numpy.diag(p**a), U.T])
    return Pa

def psi_molecule_from_file(f, frm=None, no_com=True, no_reorient=True):
    "Construct psi4.core.Molecule object from structure file"
    if frm is None: frm = f.split('.')[-1].lower()
    #
    if frm == 'xyz':
       qmol = psi4.qcdb.Molecule.init_with_xyz(f, no_com=no_com, no_reorient=no_reorient)  
       mol  = psi4.geometry(qmol.create_psi4_string_from_molecule())
    else: raise ValueError("Unrecognised format - %s -" % frm)
    #
    mol.update_geometry()
    return mol
