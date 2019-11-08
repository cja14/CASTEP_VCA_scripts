# CASTEP_VCA_scripts

Scripts for interacting with CASTEP input and output files for simulations that model solid solutions using the virtual crystal approximation (VCA).

These scripts may be used to construct input files for simulations and extract info (structures, calculation parameters, forces, energies etc) from input and output files, even for cases that use the VCA -- something that the _ase_ (atomic simulation enviroment) module cannot do. These scripts may also be used to easily convert between solid solution ("mixed") structures and "pure" structures -- which may be interpretted by standard atomistic manipulation tools. There is also the functionality (within the _phonon_VCA.py_ module) to compute harmonic phonons using _Phonopy_ for cells with mixed atoms on a single site.

The main modules, containing functions and classes for general use, are:
* _readmixcastep.py_ -- for reading/writing VCA CASTEP input/output files
* _mixmap.py_ -- for managing mapping between pure and mix structures
* _phonons_VCA.py_ -- a wrapper to manage phonon calculations with the VCA

The following modules then provide more general utilities:
* _strindices.py_ -- for identifying lines in files containing various combinations of strings
* _casase.py_ -- a wrapper to the _ase.io.read()_ method to suppress unnecessary output if CASTEP is not integrated to run within ase (e.g. if simulations are run externally).

Finally, the following scripts are command line tools for quickly manipulating structures:
* _SS_to_endmember.py_ -- converts a mixed (solid solution) structure quickly to a pure structure.
* _pinchposns.py_ -- takes two cell files as input and shifts atoms in the second to be closest to those in the first (considering PBCs).
* _shiftcas.py_ -- shifts all atoms in a cell by a given vector.

The _unittests.py_ script can be run to check that most functionality within these modules works (i.e. does not throw up an error). However for a tutorial on how to use the scripts, I recommend looking in the _examples/_ directory.