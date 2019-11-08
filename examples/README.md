# Tutorial

For this tutorial, it will be assumed that you are working in the _examples/_ directory and have the python following packages installed:

* _numpy_
* _ase_
* _phonopy_

Also, it is assumed that you have the following files in your python path:

* _readmixcastep.py_
* _mixmap.py_
* _phonons_VCA.py_

## Reading in a solid solution castep

First, load the readmixcastep module and use this to read in output form a CASTEP DFT calculation on a system using the VCA (with mixed atoms on a given site): 

```python
import readmixcastep as rc

cas = rc.readcas('Ca2.15Sr0.85Ti2O7_Amam.castep')
mixatoms = cas.extract_struc()
```



