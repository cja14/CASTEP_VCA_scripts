#!/usr/bin/env python3

import sys
import mixmap
import readmixcastep as rc
import numpy as np

cellfile = sys.argv[1]
shift = np.array([float(s) for s in sys.argv[2:5]])  # frac coords

cas = rc.readcell(cellfile)
mixatoms = cas.extract_struc()
mixkey = cas.get_mixkey()
mapping = mixmap.mixmap(mixatoms, mixkey)
posns = mixatoms.get_scaled_positions()
posns += shift
mixatoms.set_scaled_positions(posns)
mapping.casprint(mixatoms, cellfile)
