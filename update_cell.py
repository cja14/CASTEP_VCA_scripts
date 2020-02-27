#!/usr/bin/env python3

import readmixcastep as rc
import mixmap
from glob import glob

def main():

  casfiles = glob("*.castep")

  for file in casfiles:
    cellname=file.replace(".castep", ".cell")

    cas = rc.readcas(file)
    kpoints, offset = cas.get_kpoints()
    constrs=cas.get_cell_constrs()
    psps = cas.get_psps()

    mixatoms = cas.extract_struc()
    mixkey = cas.get_mixkey()

    #Map 
    mapping = mixmap.mixmap(mixatoms, mixkey)
    mapping.setcellparams(kpoints=kpoints, kpoints_offset=offset,\
        cell_constrs=constrs, pseudos=psps)

    mapping.casprint(mixatoms, cellname, pure=False)
    

if __name__ == "__main__":
  main()
