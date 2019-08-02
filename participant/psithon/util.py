#!/usr/bin/python3
"""
 Workshop Utilities
"""
import numpy
import psi4

def matrix_power(P, a):#TODO
    "Power of a real symmetric matrix P^a"
    raise NotImplementedError

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

def two_index_transform(int_ab, C1, C2):
    int_Ib = numpy.einsum("ab,aI->Ib", int_ab, C1); del int_ab
    int_IJ = numpy.einsum("Ib,bJ->IJ", int_Ib, C2); del int_Ib
    return int_IJ

def two_index_transform_full(int_ab, C1, C2):
    int_IJ = numpy.einsum("ab,aI,bJ->IJ", int_ab, C1, C2)
    return int_IJ

def four_index_transform(eri_abcd, C1, C2, C3, C4):#TODO
    raise NotImplementedError

def four_index_transform_full(eri_abcd, C1, C2, C3, C4):#TODO
    raise NotImplementedError
