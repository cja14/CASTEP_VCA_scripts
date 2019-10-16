#!/usr/bin/python

import sys
import mixmap
import readmixcastep as rc
import numpy as np

cellfile = sys.argv[1]
shift = np.array([float(s) for s in sys.argv[2:5]])  # frac coords

cas = rc.readcell(cellfile)
casatoms = cas.extract_struc()
mixkey = cas.get_mixkey()
mapping = mixmap.mixmap(casatoms, mixkey)
posns = casatoms.get_scaled_positions()
posns += shift
casatoms.set_scaled_positions(posns)
mapping.casprint(casatoms, cellfile)
