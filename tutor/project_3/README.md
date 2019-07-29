# Molerucal Dynamics

## Newtonian Molecular Dynamics

## Quantum Dynamics

## Trajectory Surface Hopping

## Algorithm

# Programming Task

Our task is to write an application that performs the TSH dynamics for a
molecular aggregate with one molecule being described by the high level
quantum mechanical method allowing for computing electronic excited states
(such as CIS, TDDFT, EOM-CC, CASSCF, ADC etc.), an the rest by a low level
ground state method (either HF, MP2, CC or just molecular mechanics). We wish to
be able to choose between various embedding schemes (electronic embedding due to constant
or variable atomic charges, electronic density, polarizable density embedding and so on).
Therefore, the code has to provide a unified and flexible framework and platform 
for extending its functionalities without disturbing the core idea and structure of the application.
It should also allow to specify the interface between third-party Quantum Chemistry software
if necessary.

## Specification of Objects and Their Relationships

The most convenient way to start is to establish the structure of the project: data types, 
objects and relationships between them. Since our task is to model the molecular dynamics, our target object
will be a *dynamical system*. The schematic representation of such a system for the need of TSH algorithm
is shown below:

<img src="../../doc/figures/plan.png" height="300"/>

The circles represent certain classes of objects understood as basic abstract building blocks
of our program.
It is important to distinguish here between the abstract classes and the actual utility classes
which will be later on used to instantiate the objects. 

The six major abstract base classes were identified:
 1. **Aggregate** - this class describes the entire molecular system. It
    contains information about the molecular composition of the aggregate, its 
    coordinates, charge and multiplicity division
    into fragments, as well as which molecule is in excited state allowed region.
 2. We probably want the instances of class Aggregate to store objects of type **Molecule**, 
    referring to particular molecules from the aggregate (such as `psi4.core.Molecule` class).
 3. **Hamiltonian** - is to describe the quantum Hamiltonian of an entire aggregate.
    That is to say, here the environmental effects on each molecular wavefunction
    are to be implemented. Object of this class will therefore be responsible for
    solving Schrodinger equation for each fragment in the presence of other fragment,
    resulting in an effective electronic Hamiltonian of an entire system.
    We probably want to refer inside its instance to Aggregate as a reference.
 4. Each Molecule object will be handled by a separate **Computer** object, that stores
    the actual information about the molecule and its wavefunction. Computer therefore can be viewed
    as sort of interface between the TSH application and certain 
    method of solving the Schrodinger equation, implemented elsewhere.
 5. Computer instances will store the information of the multireference wavefunction
    of a given fragment in instances of type **CIWavefunction**. This class
    describes the composition of an arbitrary set of electronic states in terms of Slater determinants
    as basis functions.
 6. **Trajectory** object will manage phase space elements in time. It should store 
    nuclear positions, velocities and forces, as well as copies of
    CIWavefunction objects at a given time.
    
    
## Implementation

Let us start from the end: implement the classes of smallest range from the above scheme
(Aggregate, TimePoint),
and finish with class of largest range (System).

### Aggregate

To set up high level and bath molecules, we initialize the object by constructing 
copies of each fragment and saving it in e.g. `qm` and `bath` variables.
It is also useful to define two functionalities of Aggregate instances:
 * possibility to update from an input coordinate matrix (such as `psi4.core.Matrix` or `numpy.ndarray`)
 * possibility to save current coordinates to a trajectory file

```python
class Aggregate:
  def __init__(self, psi4_molecule):
      self.all = psi4_molecule
      self.nfrags = psi4_molecule.nfragments()
      self.qm = psi4_molecule.extract_subsets(1)
      self.bath = [] if self.nfrags == 1 else [psi4_molecule.extract_subsets(2+i) for i in range(self.nfrags-1)]

  def update(self, xyz):
      "Update Aggregate coordinates"
      pass

  def save_xyz(self, out):
      "Save current coordinates to xyz trajectory file"
      pass
```

Try to figure out the implementation of the instance methods above.

### TimePoint

Implementation of this class is also relatively straightforward. 
We want to store elements of classical phase space (positions, velocities and
forces on the nucleis). We also want to store multiconfigurational wavefunction.
It is quite useful to add a functionality that saves certain data to trajectory `dat` file
(e.g., velocities, forces, electronic state, etc.).

```python
class TimePoint:
  def __init__(self, x, v, f, s):
      self.x = x.copy()
      self.v = v.copy()
      self.f = f.copy()
      self.s = s
      self.time_id = 0
  def save(self, out):
      pass
```


### Trajectory


The next class element, that is probably the smallest in range is Trajectory,
because it encompasses only TimePoint. Here we decide to store in memory
only two last time points. There is also a few of very useful functionalities
Trajectory object could have:
 * adding next point
 * setting two last points
 * saving last point to dat file
 * computing kinetic energy of last point
 * readjusting velocities (momentum rescaling)
 * canonicalization of velocities (for example, at the beginning of simulation).

```python
class Trajectory(ABC):
  def __init__(self, molecules):
      ABC.__init__(self)
      self.natoms = molecules.natom()
      self.point_last = None
      self.point_prev = None

  def add_point(self, point):
      "Add next point to the trajectory"
      pass

  def set_points(self, prev, last):
      "Set previous and current (last) point in the trajectory"
      pass

  def save(self, out): 
      "Save trajectory current state to dat file"
      self.point_last.save(out)

  def kinetic_energy(self):
      "Compute kinetic energy of nuclei"
      pass

  def rescale_velocities(self, e_kin):
      "Readjust velocities to reproduce given kinetic energy"
      pass

  def canonicalize_velocities(self, temp): 
      "Make canonical distribution of velocities for a given temperature"
      pass
```

To design the Hamiltonian class, it is best to design first Computer class, since it is of
lesser range.

### Computer

This time, we decide to create many different types of Computer, depending
on the method chosen to compute energy and forces.
Let start from the abstract base *Computer*.

```python
class Computer(ABC):
  def __init__(self, molecule):
      ABC.__init__(self)
      #
      self.molecule = molecule
      self.nuclear_repulsion_energy = molecule.nuclear_repulsion_energy()
      #
      self.ciwfn = None
      self.nstates = None
      self.forces = None

  @classmethod
  def create(cls, method, molecule, nstates=1):
      m = method.lower()
      if m == "psi4scf": return psi4SCF_Computer(molecule)
      elif m == "mycis": return myCIS_Computer(molecule, nstates)
      else: raise ValueError("Wrong method chosen for computer")

  def update(self, xyz):
      self.molecule.set_geometry(xyz)
      self.molecule.update_geometry()
      self.nuclear_repulsion_energy = self.molecule.nuclear_repulsion_energy()

  def compute(self):
      self.ciwfn = self._compute_energy()
      self.forces = self._compute_forces()

  @abstractmethod
  def _compute_energy(self): pass

  @abstractmethod
  def _compute_forces(self): pass
```

Here, we decided to implement generation of CIWavefunction and forces
in derived subclasses. For example, let us think we want to perform
ground state simulation with Hartree-Fock theory. In Psi4, HF solver
enables efficient calculation of total energy and forces.
Therefore, one of applied classes could look like this:

```python
class psi4SCF_Computer(Computer):
  def __init__(self, molecule): 
      Computer.__init__(self, molecule)
      self.nstates = 1
  def _compute_energy(self):
      e, w = psi4.energy("SCF", molecule=self.molecule, return_wfn=True)
      E = numpy.array([e - self.molecule.nuclear_repulsion_energy()])
      state = HF_CIWavefunction(w, E)
      return state
  def _compute_forces(self):
      g = -psi4.gradient("SCF", molecule=self.molecule).to_array(dense=True)
      return numpy.array([g])
```
We set number of states equal to just 1 because there is only ground state.
Note, how we compute energy and forces directly by interfacing our TSH application
with Psi4.

In order to implement some description of excited states, we choose to use our 
CIS implementation from Project II. Here, we decide to introduce several intermediate abstract classes
just to make the code even more encapsulated.
For example, we assume that so far we have no analytical forces implemented,
hence we resort to numerical calculation for all excited state methods.
This feature can be implemented in intermediate class `ExcitedState_Computer`.
Our actual applied class is realized as the last element of the following class three:

```python
class ExcitedState_Computer(Computer):
  def __init__(self, molecule, nstates):
      Computer.__init__(self, molecule)
      self.nstates = nstates

  def _compute_forces(self):
      "Numerically compute forces"
      pass

class CIS_Computer(ExcitedState_Computer):
  def __init__(self, molecule, nstates):
      ExcitedState_Computer.__init__(self, molecule, nstates) 

class myCIS_Computer(CIS_Computer):
  def __init__(self, molecule, nstates):
      CIS_Computer.__init__(self, molecule, nstates)

  def _compute_energy(self):
      cis = CIS.create(self.molecule, verbose=False, save_states=self.nstates-1, reference='uhf')
      cis.run()
      state = CIS_CIWavefunction(cis.ref_wfn, cis.E, cis.W)
      return state
```
