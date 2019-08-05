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
      self._bfs =
      # Wavefunction
      self._wfn = 
      # Number of alpha electrons
      self._ndocc = 
      # Integral calculator
      self._mints = 
      # JK object
      self._jk = 
      ### Accessors
      # nuclear repulsion energy
      self.e_nuc = 
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
      self.S = 
      self.X = self._orthogonalizer(self.S)
      return

  def run(self, maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True):#TODO
      "Solve SCF (public interface)"
      if guess is None:
         # Form Hcore                    
         ...
      else: H = numpy.asarray(guess)
      self.H = H.copy()
      self._run(H, maxit, conv, damp, ndamp, verbose)
      return

  # --- protected --- #

  def _run(self, H, maxit, conv, damp, ndamp, verbose):#TODO
      "Solve SCF (protected interface)"
      # [1] Guess density matrix
      
      # [2] Start iteration cycles

      while (abs(e_old - e_new) > conv):
        niter += 1
        # [3] form Fock matrix

        # [4] compute total energy

        if verbose:
            print (" @SCF Iter {:02} E = {:14.8f}".format(niter, e_new))

        # [5] transform Fock matrix to orthogonal AO basis           

        # [6] diagonalize the Fock matrix

        # [7] convert LCAO-MO coefficiets to non-orthogonal AO basis

        # [8] form density matrix

        # [9] save current data

        if niter > maxit: break
      return

  def _orthogonalizer(self, S):#TODO
      "Form orthogonalizer"
      return NotImplementedError
