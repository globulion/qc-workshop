# Configuration Interaction Singles

One of the easiest ways to model electronic excited states is through
the configuration interaction singles (CIS) method. 
In this method, the electronic wavefunction of a molecule
is assumed to be given by

<img src="../../doc/figures/equations/cis-ansatz.png" height="60"/>

where the first term is the Hartree-Fock (reference) Slater determinant
composed of the HF occupied molecular orbitals,
whereas the latter terms constitute of the Slater determinants
in which one occupied MO was replaced with one virtual MO, weighted by
CIS coefficients *c*. In general, any Slater determinant can be
fully described by the set of alpha and beta spin-orbitals

<img src="../../doc/figures/equations/slater-determinant-formula.png" height="40"/>

which are to be arranged in a form of a determinant
so that all the permutations of electron labels are exhausted,
and are therefore by construct antisymmetric with respect to interchange
of electron labels:

<img src="../../doc/figures/equations/slater-determinant-formula-determinant.png" height="70"/>

Effective Hamiltonian can be constructed in the basis of all these
Slater determinants. The form before integrating out the spin coordinates
is as follows:

<img src="../../doc/figures/equations/cis-hiajb.png" height="30"/>

where the Fock matrix elements are given by

<img src="../../doc/figures/equations/fock.png" height="60"/>

In the above equation, *h* is the one-electron core Hamiltonian
operator

<img src="../../doc/figures/equations/h-core.png" height="50"/>

composed of the electronic kinetic energy operator *T* and
electron-nucleus interaction energy operator *V* (the summation over *u* runs through
all nuclei in the molecule).

Note that we need to explicitly care about electron spins when evaluating
Hamiltonian matrix elements in the Slater determinant basis.
The CIS Hamiltonian can be drawn like this:

<img src="../../doc/figures/cis-hamiltonian.png" height="190"/>

The first diagonal block is just one number, equal to the reference (ground state) electronic
energy, obtained by solving the Hartree-Fock problem. The offdiagonal blocks
in white are zero due to de Bruillouin theorem. The next diagonal blocks (blue)
correspond to the excitations of electron of the same spin in bra's and ket's,
and they are given by

<img src="../../doc/figures/equations/cis-aa.png" height="40"/>

and analogously for the beta spin.
The offdiagonal blocks (green) are given by

<img src="../../doc/figures/equations/cis-ba.png" height="40"/>

and analogously for the twin term. In the above working equations,
two-electron integrals in MO basis were represented in the Coulomb (or chemist's)
notation which is given by


