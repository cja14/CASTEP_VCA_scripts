#!/usr/bin/env python3

import mixmap
import readmixcastep as rc
import sys

""" Command line tool to quickly get a pure cell file from a solid soln. """

# Load the input file (accepts .castep or .cell)
casfile = str(sys.argv[1])
if '.castep' in casfile:
    chem = casfile.replace('.castep', '')
    cas = rc.readcas(casfile)
elif '.cell' in casfile:
    chem = casfile.replace('.cell', '')
    cas = rc.readcell(casfile)

# Extract calculation parameters and convert from mix to pure structure
mixatom = cas.extract_struc(iteration=None)
mixkey = cas.get_mixkey(iteration=None)
press = cas.get_ext_press()
kpoints, offset = cas.get_kpoints()
pseudos = cas.get_psps()
constraints = cas.get_cell_constrs()
mapping = mixmap.mixmap(mixatom, mixkey)
pureatom = mapping.mix2pure(mixatom)

# Write the pure cell (always with the same handle plus the suffix _nonSS)
mapping.setcellparams(pressure=press, cell_constrs=constraints,
                      pseudos=pseudos, kpoints=kpoints, kpoints_offset=offset)
cellfile = chem+'_nonSS.cell'
mapping.casprint(pureatom, cellfile, pure=True)
