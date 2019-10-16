#!/usr/bin/python

import ase
import mixmap
import numpy as np
import readmixcastep as rc
import sys
import glob

# Script to produce a pure compound cell file from a solid solution

# WARNING: will not work with spins

casfile = str(sys.argv[1])     # First argument is .castep file
if '.castep' in casfile:
    chem = casfile.replace('.castep', '')
    cas = rc.readcas(casfile)
elif '.cell' in casfile:
    chem = casfile.replace('.cell', '')
    cas = rc.readcell(casfile)
# it = -1
it = None  # WARNING: this will look at geometries of incomplete calculations!
casatom = cas.extract_struc(iteration=it) 
mixkey = cas.get_mixkey(iteration=it)
press = cas.get_ext_press()
kpoints, offset = cas.get_kpoints()
pseudos = cas.get_psps()
constraints = cas.get_cell_constrs()
mapping = mixmap.mixmap(casatom, mixkey)
phonatom = mapping.cas2phon(casatom)
mapping.setcascellinfo(pressure=press, cell_constrs=constraints,
                       pseudos=pseudos, kpoints=kpoints, kpoints_offset=offset)
cellfile = chem+'_nonSS.cell'
mapping.casprint(phonatom, cellfile, phon=True)
