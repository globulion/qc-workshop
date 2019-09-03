#!/usr/bin/python3
"""
 Demonstrates the use of Psi4 from Python level.
 Useful notes:
  o Use psi4.core module for most of the work
  o Useful modules within psi4.core: 
     - MintsHelper
     - Molecule
     - BasisSet
     - ExternalPotential
     others
  o Psi4 defines its own matrix type (psi4.core.Matrix). 
    Extracting numpy.array is easy:
       numpy_array = numpy.asarray(psi4_matrix)
    Creating Psi4 matrix from array is also easy:
       psi4_matrix = psi4.core.Matrix.from_array(numpy_array)     
  o To compute 1-el potential matrix for a set of charges 
    use ExternalPotential (charge positions are to be provided in Angstroms)
    unless charges are just nuclei within the basis set (in this case use of ao_potential method
    of MintsHelper is easier).  
  o ao_potential method of MintsHelper is limited only for nuclei within the same basis set
    (the nuclei are taken from the first basis set axis, for example:
      mints = MintsHelper(basis_X)         
      mints.ao_potential()                 -> nuclei taken from basis of mints object (basis_X)
      mints.ao_potential(basis_1, basis_2) -> nuclei taken from basis_1
  o Psi4 has efficient and easy to use method of defining fragments within a molecule (use '--' separator).
    Defining ghost atoms and extracting fragment i in the multimer-centred basis set is also very straighforward
    (method extract_subsets(...) of psi4.core.Molecule)
---
 Bartosz Błasiak
"""

import psi4
import numpy

MAX_NBF = 128


class SCF:
  """
 ---------------------------------------------------------------------------------------------------------------
                              Self-Consistent Field (SCF) Procedure for Hartree-Fock Model
 ---------------------------------------------------------------------------------------------------------------

 Demo for RHF-SCF method (closed shells). Implements SCF algorithm
 with primitive damping of the AO Fock matrix. 

 Usage:
  scf = SCF(molecule)
  scf.run(maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True)

 The above example runs SCF on 'molecule' psi4.core.Molecule object 
 starting from core Hamiltonian as guess (guess=None) 
 and convergence 1.0E-7 A.U. in total energy with 30 maximum iterations
 (10 of which are performed by damping of the Fock matrix with damping coefficient of 0.01).
 The SCF iterations are printed to standard output (verbose=True).
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, May 4th 2018
"""
  def __init__(self, mol):
      "Initialize BasisSet, Wavefunction and JK objects"
      # Basis set
      self._bfs = psi4.core.BasisSet.build(mol, "ORBITAL", psi4.core.get_global_option("BASIS"), 
                                                  puream=psi4.core.get_global_option("PUREAM"))
      # Number of alpha electrons
      q_el = int(mol.molecular_charge() - sum([mol.Z(x) for x in range(mol.natom())]))
      if q_el%2 == 1: raise ValueError("The molecule has an unpaired electron! Only closed-shells are supported")
      self._ndocc =int(-q_el / 2)
      # Integral calculator
      self._mints = psi4.core.MintsHelper(self._bfs)
      # JK object
      self._jk = psi4.core.JK.build(self._bfs, jk_type="Direct")
      self._jk.set_memory(int(5e8))
      self._jk.initialize()
      ### Accessors
      # nuclear repulsion energy
      self.e_nuc = mol.nuclear_repulsion_energy()
      # Total Energy
      self.E = None
      # Density Matrix
      self.D = None
      # LCAO-MO coeffs (occ)
      self.Co= None
      # LCAO-MO coeffs (occ+vir)
      self.C = None
      # Fock matrix 
      self.F = None
      # Orbital energies
      self.eps = None
      # Hcore matrix
      self.H = None
      # Overlap integrals and orthogonalizer
      self.S = numpy.asarray( self._mints.ao_overlap() )
      self.X = self._orthogonalizer(self.S)
      return

  def run(self, maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, v_ext=None):
      "Solve SCF (public interface)"
      # Form Hcore                    
      T = self._mints.ao_kinetic()
      V = self._mints.ao_potential()
      H = T.clone()
      H.add(V)
      H = numpy.asarray(H)
      if v_ext is not None:
         H += numpy.asarray(v_ext)
      self.H = H.copy()
      # Getermine guess
      if guess is None:
         guess = H
      else:
         guess = numpy.asarray(guess)
      # Run SCF
      self._run(guess, maxit, conv, damp, ndamp, verbose)
      return

  # --- protected --- #

  def _run(self, guess, maxit, conv, damp, ndamp, verbose):
      "Solve SCF (protected interface)"
      # First step: Guess density matrix
      F    = numpy.dot(self.X, numpy.dot(guess, self.X)) 
      E, C = numpy.linalg.eigh(F)
      C    = numpy.dot(self.X, C)
      idx = numpy.argsort(E)
      E = E[  idx]
      C = C[:,idx]
      Co= C[:,:self._ndocc].copy()
      D = numpy.dot(Co,Co.T)
      
      niter = 0
      e_old = 1e8
      e_new = 1e7
      F_old = guess.copy()
      while (abs(e_old - e_new) > conv):
        niter += 1
        # form Fock matrix
        self._jk.C_clear()
        self._jk.C_left_add(psi4.core.Matrix.from_array(Co, "C matrix"))
        self._jk.compute()
        F_new = self.H + 2.0 * numpy.asarray(self._jk.J()[0]) - numpy.asarray(self._jk.K()[0])
        if niter < ndamp: 
           F = damp * F_old + (1.0 - damp) * F_new
        else:             
           F = F_new
        F_old = F.copy()
        # compute total energy
        e_old = e_new
        e_new = numpy.trace( numpy.dot(D, self.H + F) ) + self.e_nuc
        if verbose:
            print (" @SCF Iter {:02} E = {:14.8f}".format(niter, e_new))
        # transform Fock matrix to orthogonal AO basis           
        F = numpy.dot(self.X, numpy.dot(F, self.X))
        # diagonalize the Fock matrix
        E, C = numpy.linalg.eigh(F)
        # convert LCAO-MO coefficiets to non-orthogonal AO basis
        C = numpy.dot(self.X, C)
        # form density matrix
        idx = numpy.argsort(E) 
        E = E[  idx]
        C = C[:,idx]
        Co= C[:,:self._ndocc].copy()
        D = numpy.dot(Co,Co.T)
        # save
        self.D = D.copy()
        self.E = e_new
        self.F = F_old.copy()
        self.C = C.copy()
        self.Co= Co.copy()
        self.eps = E.copy()
        if niter > maxit: break
      return

  def _orthogonalizer(self, S):
      "Form orthogonalizer"
      L, U = numpy.linalg.eigh(S)
      L    = numpy.diag(1./numpy.sqrt(L))
      X    = numpy.dot(U, numpy.dot(L, U.T))
      return X

