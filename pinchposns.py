#!/usr/bin/python

import sys
import mixmap
import readmixcastep as rc
import numpy as np

"""
Reads two cell files and ensures that the positions in the second are close
to those in the first.
"""

cellfile0 = sys.argv[1]
cellfile1 = sys.argv[2]

###########################################################################


def reorder_atoms(LSatoms, HSatoms):
    N = len(HSatoms)
    if len(LSatoms) != N:
        raise ValueError('Low symmetry structure has different number of ' +
                         'atoms to high-symmetry strucutre.')
    LSposns = LSatoms.get_positions()
    HSposns = HSatoms.get_positions()
    HSelems = HSatoms.get_chemical_symbols()
    cell = LSatoms.get_cell()
    # For the moment I won't match elems
    LSposnew = np.zeros((N, 3))
    for i, HSpos in enumerate(HSposns):
        ds = np.zeros((N))
        for j, jpos in enumerate(LSposns):
            ds[j] = mixmap.closestdist(HSpos, jpos, cell)
        j = np.where(ds == min(ds))[0][0]
        LSposnew[i, :] = mixmap.closestimage(HSpos, LSposns[j], cell)
    LSatoms.set_chemical_symbols(HSelems)
    LSatoms.set_positions(LSposnew)
    return LSatoms


##########################################################################

if '.cell' in cellfile0:
    cas0 = rc.readcell(cellfile0)
elif '.castep' in cellfile0:
    cas0 = rc.readcas(cellfile0)
casatoms0 = cas0.extract_struc()
mixkey0 = cas0.get_mixkey()
mapping0 = mixmap.mixmap(casatoms0, mixkey0)
phonatoms0 = mapping0.cas2phon(casatoms0)
posns0 = phonatoms0.get_positions()
Nions = len(phonatoms0)

cas1 = rc.readcell(cellfile1)
casatoms1 = cas1.extract_struc()
mixkey1 = cas1.get_mixkey()
mapping1 = mixmap.mixmap(casatoms1, mixkey1)
phonatoms1 = mapping1.cas2phon(casatoms1)
posns1 = phonatoms1.get_positions()
cell1 = phonatoms1.get_cell()

spins = cas1.get_init_spin()
kpoints, offset = cas1.get_kpoints()
constrs = cas1.get_cell_constrs()
psps = cas1.get_psps()
pressure = cas1.get_ext_press()

if len(phonatoms1) != Nions:
    raise ValueError(cellfile0 + ' has ' + str(Nions) +
                     ' distinct atoms whereas ' + cellfile1 + ' has only ' +
                     str(len(phonatoms1)))

phonatoms1 = reorder_atoms(phonatoms1, phonatoms0)

casatoms1 = mapping1.phon2cas(phonatoms1)
mapping1.setcascellinfo(pseudos=psps, kpoints=kpoints, kpoints_offset=offset,
                        spins=spins, pressure=pressure, cell_constrs=constrs)
mapping1.casprint(casatoms1, cellfile1)

