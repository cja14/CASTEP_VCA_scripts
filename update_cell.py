#!/usr/bin/env python3
"""
This script converts a set of geometry optimised .castep output files of solid
solutions to pure .cell input files.
"""

import readmixcastep as rc
import mixmap
from glob import glob
from convert2cell import cell2cif

def main(pureelems):

    casfiles = glob("*.castep")

    for cfile in casfiles:
        print(cfile)
        cellname=cfile.replace(".castep", ".cell")
        if not glob(cellname):
            cas = rc.readcas(cfile)
            kpoints, offset = cas.get_kpoints()
            constrs=cas.get_cell_constrs()
            psps = cas.get_psps()

            mixatoms = cas.extract_struc()
            mixkey = cas.get_mixkey(pureelems=pureelems)

            #Map 
            mapping = mixmap.mixmap(mixatoms, mixkey)
            pureatoms = mapping.mix2pure(mixatoms)

            mapping.casprint(pureatoms, cellname, pure=True)
            cell2cif(cellname)
        elif not glob(cellname.replace(".cell", ".cif")):
            cell2cif(cellname)

def mainpure():
    casfiles = glob("*.castep")
    for cfile in casfiles:
        print(cfile)
        cellname=cfile.replace(".castep", ".cell")
        if not glob(cellname):
            cas = rc.readcas(cfile)
            kpoints, offset = cas.get_kpoints()
            constrs=cas.get_cell_constrs()
            psps = cas.get_psps()

            pureatoms = cas.extract_struc()
            mixkey = cas.get_mixkey()

            mapping = mixmap.mixmap(pureatoms, mixkey)
            mapping.setcellparams(pseudos=psps, kpoints=kpoints,\
                kpoints_offset=offset, cell_constrs=constrs)
            mapping.casprint(pureatoms, cellname, pure=True)
        elif not glob(cellname.replace(".cell", ".cif")):
            cell2cif(cellname)

    return None

def mix(purelems):

    casfiles = glob("*.castep")

    for cfile in casfiles:
        print(cfile)
        cellname=cfile.replace(".castep", ".cell")
        if not glob(cellname):
            cas = rc.readcas(cfile)
            kpoints, offset = cas.get_kpoints()
            constrs=cas.get_cell_constrs()
            psps = cas.get_psps()

            mixatoms = cas.extract_struc()
            mixkey = cas.get_mixkey(pureelems=pureelems)

            #Map 
            mapping = mixmap.mixmap(mixatoms, mixkey)

            mapping.casprint(mixatoms, cellname, pure=False)


if __name__ == "__main__":
    #pureelems = {'Sr': 'La', 'Al': 'Mg', 'Ca': 'La', 'Ba': 'La', 'Nd': 'La',\
    #            'Y': 'La'}
    mainpure()
