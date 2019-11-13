#!/usr/bin/env python3

import sys
import mixmap
import readmixcastep as rc
import numpy as np

"""
Reads two cell files and ensures that the positions in the second are close
to those in the first.
"""


def reorder_atoms(fixatoms, pinchatoms):
    """ Reorder pinchatoms (i.e. accounting for PBCs) to match fixatoms
    Note: this only works for pure atomic structures (no mixed atoms) """
    N = len(fixatoms)
    if len(pinchatoms) != N:
        raise ValueError('Low symmetry structure has different number of ' +
                         'atoms to high-symmetry strucutre.')
    pinchposns = pinchatoms.get_positions()
    fixposns = fixatoms.get_positions()
    fixelems = fixatoms.get_chemical_symbols()
    cell = pinchatoms.get_cell()
    # For the moment I won't match elems
    pinchposnew = np.zeros((N, 3))
    for i, fixpos in enumerate(fixposns):
        ds = np.zeros((N))
        for j, jpos in enumerate(pinchposns):
            bprime, dist = mixmap.mixmap.closestimage(fixpos, jpos, cell,
                                                      rtn_dist=True)
            ds[j] = dist
        j = np.where(ds == min(ds))[0][0]
        pinchposnew[i, :] = mixmap.mixmap.closestimage(fixpos, pinchposns[j],
                                                       cell, rtn_dist=False)
    pinchatoms.set_chemical_symbols(fixelems)
    pinchatoms.set_positions(pinchposnew)
    return pinchatoms


##########################################################################

if __name__ == '__main__':
    """ Run from the command line, pinch the second cell to match the first """
    cellfile0 = sys.argv[1]
    cellfile1 = sys.argv[2]

    # Load the fixed cell (incl generating pure structure)
    if '.cell' in cellfile0:
        mix0 = rc.readcell(cellfile0)
    elif '.castep' in cellfile0:
        mix0 = rc.readcas(cellfile0)
    mixatoms0 = mix0.extract_struc()
    mixkey0 = mix0.get_mixkey()
    mapping0 = mixmap.mixmap(mixatoms0, mixkey0)
    pureatoms0 = mapping0.mix2pure(mixatoms0)
    posns0 = pureatoms0.get_positions()
    Nions = len(pureatoms0)
    
    # Load the cell to be pinched (incl generating pure structure)
    mix1 = rc.readcell(cellfile1)
    mixatoms1 = mix1.extract_struc()
    mixkey1 = mix1.get_mixkey()
    mapping1 = mixmap.mixmap(mixatoms1, mixkey1)
    pureatoms1 = mapping1.mix2pure(mixatoms1)
    posns1 = pureatoms1.get_positions()
    cell1 = pureatoms1.get_cell()
    
    # Load calculation params for the second cell (for writing new cell file)
    spins = mix1.get_init_spin()
    kpoints, offset = mix1.get_kpoints()
    constrs = mix1.get_cell_constrs()
    psps = mix1.get_psps()
    pressure = mix1.get_ext_press()
    
    # Reorder the positions of the pure structure
    pureatoms1 = reorder_atoms(pureatoms0, pureatoms1)
    
    # Transform back to the mixed structure and write the cell file
    mixatoms1 = mapping1.pure2mix(pureatoms1)
    mapping1.setcellparams(pseudos=psps, kpoints=kpoints,
                           kpoints_offset=offset, spins=spins,
                           pressure=pressure, cell_constrs=constrs)
    mapping1.casprint(mixatoms1, cellfile1)
