# Tutorial

For this tutorial, it will be assumed that you are working in the _examples/_ directory and have the python following packages installed:

* _numpy_
* _ase_
* _phonopy_

Also, it is assumed that you have the following files in your python path:

* _readmixcastep.py_
* _mixmap.py_
* _phonons_VCA.py_

## Reading in a solid solution .castep file

First, load the readmixcastep module and use this to read in output form a CASTEP DFT calculation on a system using the VCA (with mixed atoms on a given site): 

```python
import readmixcastep as rc

# read .castep file
cas = rc.readcas('Ca2.15Sr0.85Ti2O7_Amam.castep')

# Extract calculation parameters
task = cas.get_task()
kpoints, offset = cas.get_kpoints()
constrs = cas.get_cell_constrs()

# Extract physical parameters
energy = cas.get_energy()
forces = cas.get_forces(iteration=-1)
stresses = cas.get_stresses(iteration=0)
posns = cas.get_posns(iteration=None)

# Extract unit cell and site mixing weights
mixatoms = cas.extract_struc()
mixkey = cas.get_mixkey()
```

This creates a _readcas_ instance that we have called "cas". Using this object, we can extract parameters about the DFT calculation; such as the "task" (in this case a "geometry optimisation"), the k-point grid used (and the offset vector) or the cell constraints employed. We can also extract physical parameters that vary throughout the simulation, such as the cell energy, atomic forces, cell stresses or atomic positions. Note that we can use the optional "iteration" variable to specify which iteration to take data from (for a geometry optimisation for example). In this example we take energy and forces from the final (default) iteration -- we just say this explicitly with the forces, and stresses from the first iteration. If an integer value is given for iteration, then data will only be extracted from complete iterations (where the SCF cycle has converged). However, if a _None_ value is parsed, as in the example to extract positions, data will be extracted from the last iteration in the file, whether or not it completed successfully (this can be useful for example to restart a simulation that timed out).

The method _extract_struc()_ loads an _ase.Atoms_ object (here we've called this _mixatoms_) which contains all atoms in the simulation (including those appearing on the same site) as separate full atoms. You may have noticed that multiple entries appeared in the atomic forces and atomic positions, even for ions occupying the same site. This may be problematic because it creates a unit cell with more physical atoms than should be really be present in a single unit cell of the crystal.

Since many atomistic and crystallographic tools do not consider mixed atoms on a single site, for many operations it would be useful to deal with a "pure" structure with only one atom per site. We also don't want to throw away information about the mixing (as this _mixatoms_ object has done).

To account for this weighting, we also extract a _mixkey_ object from our castep file.

Inspect this _mixkey_ :

```python
print(mixkey)
{'0.752216 0.750000 0.900520': ('Ti', {'Ti': 1.0}),
'0.805341 0.250000 0.697388': ('O', {'O': 1.0}),
'0.736620 0.250000 0.814495': ('Ca', {'Ca': 0.7166, 'Sr': 0.2834}), ...
```

You see that it is a dictionary mapping each atomic site (given in fractional cell coordinates) to a tuple. The first element of this tuple is a single element (that we will later use for our "pure" structure), whereas the second element is a dictionary of one or more elements, where each element is paired with its weight for that site.

```python
import mixmap

mapping = mixmap.mixmap(mixatoms, mixkey)
pureatoms = mapping.mix2pure(mixatoms)
pureforces = mapping.mix2pure_forces(forces)
```

Combining these _mixkey_ and _mixatoms_ objects, we can can create a _mixmap_ class instance (from the module of the same name) which we will call _mapping_ . As the name suggests, once initialised, this _mapping_ object can be used to "map" between "mix" and "pure" structures. It may also be used to transform forces on all VCA atoms into a single force vector for each site, or to compute the masses of "pure" single-site atoms.

```python
mapping.setcellparams(kpoints=kpoints, kpoints_offset=offset, cell_constrs=constrs)
mapping.casprint(mixatoms, 'Ca2.15Sr0.85Ti2O7_Amam_mix.cell', pure=False)  # write mix cell
mapping.casprint(mixatoms, 'Ca2.15Sr0.85Ti2O7_Amam_pure.cell', pure=True)  # write pure cell
```

Using the _mapping_ object, we can output .cell files for our "mix" and "pure" structures using the _casprint_ method. Note that we set calculation parameters using the _setcellparams_ method.

## Generating a solid solution .cell file

The _readmixcastep_ and _mixmap_ modules can also be used the other way, to create "mixed" cells for a VCA calculation from "pure" cell structures.

For example, load a pure Ca3Ti2O7 structure again using readmixcastep (it works for "pure" compositions too) and define your own mixkey in the python script:

```python
cas = rc.readcell('Ca3Ti2O7_P42mnm.cell')
pureatoms = cas.extract_struc()

mixkey = {'Ca': {'Ca': 0.5, 'Sr': 0.5},
          'Ti': {'Ti': 1.0},
          'O': {'O': 1.0}}
```

Previously we saw a mixkey in the "site" format. We may also define a mixkey in the "element" format, where each element in the pure structure may be replaced by a solid solution. Note that this format for the mixkey would not work if several sites containing the same "pure" element map to different mixtures.


```python
mixatoms = mixmap.create_mixture(pureatoms, mixkey)
mapping = mixmap.mixmap(mixatoms, mixkey)

mapping.casprint(mixatoms, 'Ca1.5Sr1.5Ti2O7_P42mnm_mix.cell')
```

Using the _create_mixture()_ function in the _mixmap_ module (note not a method of the _mixmap_ class) we combine this pure structure and mixkey to give us a mixed structure. We may thereby create a _mapping_ object and use this mapping to write a mixed cell file.

## Phonon Calculations

Unfortunately CASTEP cannot perform phonon calculations using DFPT when the VCA is employed, therefore, to compute lattice dynamics we will have to use a finite differences method. The _phonons_VCA_ module acts as a wrapper to the _phonopy_ module, converting VCA mix structures and forces into pure structures to be analysed by _phonopy_ .

```python
import phonons_VCA as pVCA

pVCA.gen_perturbations('Ca1.5Sr0.5GeO4/Ca1.5Sr0.5GeO4.castep')
```

To generate the symmetry-independent displacements for a Gamma-point phonon calculation using this module requires only a single line of code. You should now discover that 21 .cell files, named _Ca1.5Sr0.5GeO4_*.cell_ , have now been generated in the sub-directory _Ca1.5Sr0.5GeO4/_ .

In normal use you will need to take these cell files and perform singlepoint CASTEP simulations. However, in true Blue Peter style, here this has been done for you and in this directory we can find appropriately named .castep output files.

Reading the files, extracting the forces, constructing and then diagonalising the dynamical matrix has been reduced using the _phonon_VCA_ module to a single function. There are two methods of running this function. Method 1 effectively generates a new set of displacements for the parent structure and assumes that these are labelled exactly the same as the output files being read. This method is quick but would fail if the perturbations change order. Method 2 checks what .castep files are available that follow the appropriate naming convention, reads these files and extracts the displacements. This method is thus slower but in theory more stable. In this example we use the latter method, with the "verbose" flag set to _True_ so we can monitor the progress of the calculation:

```python
output = pVCA.calc_phonons(Ca1.5Sr0.5GeO4/Ca1.5Sr0.5GeO4.castep',
                           method=2, verbose=True)

q_points, distances, frequencies, eigenvectors = output
freqs = frequencies[0][0]
eigvecs = eigenvectors[0][0]

print("\nFirst 10 frequencies:")
print(freqs[:10])

print("\nFirst 5 displacments of first two eigenvectors:")
print(eigvecs[1, :15].reshape(5, 3))
print(eigvecs[2, :15].reshape(5, 3))
```

The output to the _pVCA.calc_phonons()_ wrapper is the same as the underlying _phonopy_ method: _PhononObj.get_band_structure()_ . A tuple of four lists are returned where each list element corresponds to a phonon calculation performed at a different q-vector. These lists represent the q-vector, distance (for plotting a band structure), harmonic frequencies and phonon eigenvectors for each of these calculations. In this example, we compute only Gamma-point modes and therefore we show how to access the first 10 of the 3*N frequencies (in cm-1) and complex eigenvector components for the first 5 atoms of the first two eigenvectors out of the 3*N total 3*N-dimensional eigenvectors (if N is the number of ions in the pure structure).

