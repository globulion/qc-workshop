*****
![alt text](./doc/figures/toc.jpg "Logo Title Text 1")
*****

Quantum Chemistry Workshop
==========================

Quantum Chemistry with Python - Developing Own Scientific Ideas

 - Author: Bartosz Błasiak (blasiak.bartosz@gmail.com)
 - Venue: Wrocław University of Science and Technology, Wrocław, Poland
 - Date: 23-25 September 2019 (classroom and exact time will be set after 18th September)

## Table of Contents:
 * [Psithon: Python and Psi4 Combined](./tutor/psithon/README.md#using-psi4-as-python-module)
 * [Project I: Population Analyses](./tutor/project_1/README.md#population-analyses)
 * [Project II Hartree-Fock Model](./tutor/project_2/README.md#the-hartree-fock-model)
 * [Project III: Configuration Interaction with Singles](./tutor/project_3/README.md#configuration-interaction-singles)
 * [Project IV: Trajectory Surface Hopping Dynamics](./tutor/project_4/README.md#trajectory-surface-hopping)

## Description

This workshop is dedicated to theoretical quantym chemistry
from its very technical side: mathematical models and their
implementation in an object-oriented computer code that is
easier to maintain and more compatible with modern libraries. 

Usually, an object-oriented 
platform for developing own quantum chemistry models or codes needs to be
set up by a scientist alone, with an aid of a complicated (though extremely optimized
in terms of efficiency) quantum chemistry codes written mostly in Fortran.
Although it is absolutely crucial for a developer to understand
the construction and implementation of at least one of these packages 
(Gaussian, GAMESS, etc), it is not at all straighforward to adapt such
enormous and usually subroutine-based codes to a custom model development.

Object-oriented computer languages such as C++ offer more flexibility
in terms of code readability, interface with the user and maintainability,
and many new quantum chemistry codes are already being developed in this fashion
(Q-Chem, Psi4, etc). 

Python is a high-level (interpreter-based) scripting programming language
that offers a very powerful object-oriented functionalities
and interface with highly efficient mathematical libraries.
Moreover,
Psi4 provides an impressive interface with Python
which enables everybody to access quite involved quantum chemistry routines
directly from a Python script.
Therefore, it is a perfect tool for a quantum chemistry developers as 
well as students who want to broaden their knowledge in electronic structure theory.

After completing this workshop you will be more familiar with 
possibilities of writing quite sophisticated and powerful
Python applications that can do more than you could expect! During this 1-week 
interactive seminar session, we will focus on Python-Psi4 interface
to write a few quantum chemistry codes and simply have fun.

## Prerequisites
 * Basic knowledge of electronic structure theory and molecular dynamics
   - Schrodinger Equation
   - Hartree-Fock theory
   - Configuration Interaction theory
   - Classical Molecular Dynamics model
 * Unix or Mac operating system (for Windows: Microsoft Ubuntu)
 * [Basic knowledge of Python 3](https://www.programiz.com/python-programming/tutorial)
 * [NumPy module](https://numpy.org)
 * [Installation of Psi4](./doc/misc/psi4install.md)

## How to prepare

Revise the following material prior to the workshop:
 * Molecular orbitals and Slater determinant, LCAO-MO expansion
 * Hartree-Fock Theory (e.g., Ref. 1, Chapter III)
 * Configuration Interaction Singles (e.g., Ref. 2)
 * One- and two-electron integrals
 * Excited electronic states, Jablonski diagram, conical intersections, non-adiabatic coupling

Review the following literature: [5-7]

## References:

[1] [A. Szabo, N. S. Ostlund: Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory](https://store.doverpublications.com/0486691861.html)

[2] [David C. Sherril's Lecture Notes](http://vergil.chemistry.gatech.edu/notes/cis/cis.html)

[3] [Crawford's group Psithon Tutorial](https://github.com/CrawfordGroup/ProgrammingProjects)

[5] [Smith, D. A. G., Burns, L. A., et al.; *J. Chem. Theory Comput.* **2018** 14, 3504](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00286)

[5] [Hammes-Schiffer, S., Tully, J. C.; *J. Chem. Phys.* **1994** 101, 4657](https://aip.scitation.org/doi/10.1063/1.467455)

[6] [Granucci, G., Persico, M.; *J. Chem. Phys.* **2007** 126, 134114](https://aip.scitation.org/doi/10.1063/1.2715585)

[7] [Plasser, F., Ruckenbauer, M., Mai, S., Oppel, M., Marquetand, P., González, L.; *J. Chem. Theory Comput.* **2016** 12, 1207](https://pubs.acs.org/doi/10.1021/acs.jctc.5b01148)
