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

def two_index_transform(int_ab, C1, C2):
    int_Ib = numpy.einsum("ab,aI->Ib", int_ab, C1); del int_ab
    int_IJ = numpy.einsum("Ib,bJ->IJ", int_Ib, C2); del int_Ib
    return int_IJ

def two_index_transform_full(int_ab, C1, C2):
    int_IJ = numpy.einsum("ab,aI,bJ->IJ", int_ab, C1, C2)
    return int_IJ

def four_index_transform(eri_abcd, C1, C2, C3, C4):
    eri_Ibcd = numpy.einsum("abcd,aI->Ibcd", eri_abcd, C1); del eri_abcd  # cost: n^4 o
    eri_IJcd = numpy.einsum("Ibcd,bJ->IJcd", eri_Ibcd, C2); del eri_Ibcd  # cost: n^3 o^2
    eri_IJKd = numpy.einsum("IJcd,cK->IJKd", eri_IJcd, C3); del eri_IJcd  # cost: n^2 o^3
    eri_IJKL = numpy.einsum("IJKd,dL->IJKL", eri_IJKd, C4); del eri_IJKd  # cost: n   o^4
    return eri_IJKL

def four_index_transform_full(eri_abcd, C1, C2, C3, C4):
    eri_IJKL = numpy.einsum("abcd,aI,bJ,cK,dL->IJKL", eri_abcd, C1, C2, C3, C4) # cost: n^4 o^4
    return eri_IJKL
   
def _reorder(P,sim,axis=0):
    """Reorders the tensor according to <axis> (default is 0). 
<sim> is the list of pairs from 'order' function. 
In normal numbers (starting from 1...).
Copied from LIBBBG code."""
    P_new = numpy.zeros(P.shape,dtype=numpy.float64)
    if   axis==0:
         for i,j in sim:
             P_new[i-1] = P[j-1]
    elif axis==1:
         for i,j in sim:
             P_new[:,i-1] = P[:,j-1]
    elif axis==2:
         for i,j in sim:
             P_new[:,:,i-1] = P[:,:,j-1]
    return P_new

def dot_sim(a,b):
    return numpy.dot(a.ravel(),b.ravel()) / numpy.linalg.norm(a) / numpy.linalg.norm(b)
def dis_sim(a,b):
    return numpy.sum((a-b)**2)
def cos_sim(a,b):
    return dot_sim(numpy.cos(a), numpy.cos(b))
def _order(R,P,start=0,lprint=1):
    """order list: adapted from LIBBBG code"""
    new_P = P.copy()
    sim   = []
    rad =  []
    for i in range(len(R)-start):
        J = 0+start
        r = 1.0E+100
        rads = []
        a = R[i+start]
        for j in range(len(P)-start):
            #r_ = numpy.sum(( R[i+start]-P[j+start])**2)
            #r__= numpy.sum((-R[i+start]-P[j+start])**2)
            b = P[j+start]
            r_ = dis_sim(a,b)
            r__= dis_sim(-a,b)
            if r__<r_: r_=r__
            rads.append(r_)
            if r_<r:
               r=r_
               J = j
        sim.append((i+1,J+1))
        new_P[i+start] = P[J+start]
        rad.append(rads)
    for i in range(len(R)-start):
        s = numpy.sum(numpy.sign(new_P[i])/numpy.sign(R[i]))
        if lprint: print("%10d %f" %(i+1,s))
        r_ = sum(( R[i+start]-new_P[i+start])**2)
        r__= sum((-R[i+start]-new_P[i+start])**2)
       
        #if s < -154: 
        #   print "TUTAJ s < -154"
        #   #new_P[i]*=-1.
        if r__<r_:
          if lprint: print("    HERE r__ < r_ (sign reversal)")
          new_P[i]*=-1.
    return new_P, sim#, array(rad,dtype=float)


def rearrange_eigenpairs(u, u_ref, n=None, return_sim=False):
    """
 Rearrange eigenpairs. Also rephase eigenvectors if sign difference is detected.
 Inputs     : n: eigenvalues, u: eivenvectors (by column), u_ref: reference eigenvectors (by column)
 Requirement: u and u_ref need to be composed of same eigenvectors (sign arbitrary) that are in different order
              n and u need to have the same order of eigenelements.
 Returns    : n_new, u_new - when n is provided 
              u_new        - when n is not provided
 (optional) :
              sim          - similarity assignment list if return_sim=True. Returned as the last element.
"""
    u_new, sim = _order(u_ref.T, u.T, lprint=0)
    u_new  = u_new.T
    if n is not None:
       n_new  = _reorder(n, sim)
       if return_sim: return n_new, u_new, sim
       else: return n_new, u_new
    else:
       if return_sim: return u_new, sim
       else: return u_new

def check_sim(l):
    """check the sim list"""
    log = ' --- OK ---'
    for x,y in l:
        i=0;j=0
        for a,b in l:
            if a==x: i+=1
            if b==y: j+=1
        if (i>1 or j>1): 
            log = " --- !ERROR! --- "
            break
    return log
