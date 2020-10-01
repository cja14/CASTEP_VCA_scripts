#!/usr/bin/env python3

import readmixcastep as rc
import mixmap
from glob import glob
import re
import numpy as np

def get_structure(cellname):

  cas = rc.readcell(cellname)
  pureatoms = cas.extract_struc()
  kpoints, offset = cas.get_kpoints()
  constrs = cas.get_cell_constrs()
  psps = cas.get_psps()

  return pureatoms, kpoints, offset, constrs, psps

def generate_mix(x, mixname, pureatoms, mixkey, kpoints, offset, constrs, psps):

  mixatoms = mixmap.create_mixture(pureatoms, mixkey)
  mapping = mixmap.mixmap(mixatoms, mixkey)
  mapping.setcellparams(kpoints=kpoints, kpoints_offset=offset,\
      cell_constrs=constrs, pseudos=psps)

  print("Mixatoms: ", mixatoms)
  print("mixkey: ", mixkey)
  mapping.casprint(mixatoms, mixname)


def gen_mixkey(x, A_site_doping="Ba"):
  return {'La': {'La': round((2-x)/2, 4), A_site_doping: round(x/2, 4)},\
        'Mg': {'Mg': round(1-x, 4), 'Al': round(x, 4)},\
        'O': {'O': 1.0}}


def main():

    files = glob("*.cell")
    print(files)
    doping = [0.125]
    for cellname in files:
        pureatoms, kpoints, offset, constrs, psps = get_structure(cellname)
        print("Pureatoms:", pureatoms)
        for x in doping:
            x = round(x, 3)
            mixkey = gen_mixkey(x)
            generate_mix(x, cellname, pureatoms, mixkey, kpoints, offset, constrs, psps)

if __name__ == "__main__":
    main()


