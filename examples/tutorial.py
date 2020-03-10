#!/usr/bin/env python3

import readmixcastep as rc
from castep_bands import parse_seed
import mixmap
from ase.io import write

def main(pure=True):
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
    mapping.setcellparams(kpoints=kpoints, kpoints_offset=offset, cell_constrs=constrs)
    if pure:
        pureatoms=mapping.mix2pure(mixatoms)
        mapping.casprint(pureatoms, casseed + "_tutopure.cell", pure=pure)
    else:
        mapping.casprint()


def writeCIF():

    casfile=parse_seed()
    cifFile = casfile.replace(".castep", "_pure.cif")

    cas = rc.readcas(casfile)
    mixatoms = cas.extract_struc()
    mixkey = cas.get_mixkey()

    mapping = mixmap.mixmap(mixatoms, mixkey)
    pureatoms = mapping.mix2pure(mixatoms)
    write(cifFile, pureatoms)

    return pureatoms


    
if __name__ == "__main__":
    writeCIF()
