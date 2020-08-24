#!/usr/bin/env python3

import readmixcastep as rc
import mixmap
from glob import glob

def main(pureelems):

    casfiles = glob("*.castep")

    for file in casfiles:
        print(file)
        cellname=file.replace(".castep", ".cell")
        if not glob(cellname):
            cas = rc.readcas(file)
            kpoints, offset = cas.get_kpoints()
            constrs=cas.get_cell_constrs()
            psps = cas.get_psps()

            mixatoms = cas.extract_struc()
            mixkey = cas.get_mixkey(pureelems=pureelems)

            #Map 
            mapping = mixmap.mixmap(mixatoms, mixkey)
            pureatoms = mapping.mix2pure(mixatoms)

            mapping.casprint(pureatoms, cellname, pure=True)

if __name__ == "__main__":
    pureelems = {'Sr': 'La', 'Al': 'Mg'}
    main(pureelems)
