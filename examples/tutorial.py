#!/usr/bin/env python3

import readmixcastep as rc
from castep_bands import parse_seed
import mixmap

def main():
    #Parse .castep filename
    casfile = parse_seed()
    casseed = casfile.replace(".castep", "")

    #Read .castep file
    cas = rc.readcas(casfile)

    print("Task: ", cas.get_task())

    kpoints, offset = cas.get_kpoints()
    constrs = cas.get_cell_constrs()
    print("k-points: ", kpoints)
    print("offset: ", offset)
    print("Constraints: ", cas.get_cell_constrs())

    print("---------------------------------------")
    print("Energy: ", cas.get_energy(), " eV")


    print("---------------------------------------")
    print("Now creating the .cell file (mixed)")

    mixatoms = cas.extract_struc()
    mixkey = cas.get_mixkey()


    mapping = mixmap.mixmap(mixatoms, mixkey)
    mapping.setcellparams(kpoints=kpoints, kpoints_offset=offset,
            cell_constrs=constrs)
    mapping.casprint(mixatoms, casseed + "_tuto.cell", pure=False)

    
if __name__ == "__main__":
    main()
