"""
Module: Miscellaneous tools (misc-tools)

This module contains a set of miscellaneous tools for pre-processing input and
output files of electronic structure calculations.

Author: Christopher Keegan
Date: 25/08/2020

Functions taken and adapted from Chris. J. Ablitt: strindex, strindices

Strindex & Strindices:
Functions to detect the list indices in which strings occur.
Generally the list is normally to be a list of file lines from f.readlines().

"""

def pzero(x):
    """
    Function ensures a value of zero is displayed positive.

    Parameters:
    -----------
    x: float
        Value to set to 0.0 if zero.

    Returns:
    --------
    x: float
        Positive zero value.
    """
    if x == 0.0:
        x = 0.0
    return x

def posstring(posn):
    """
    This function unambiguously flattens a position array into a string.
    """
    return ' '.join([str('{0:.6f}'.format(pzero(posn[j])))
        for j in range(len(posn))])

def strindex(flist, strings, nmin=0, nmax=None, first=False, either=False):
    """ Extract index of first/last line in list at which string occurs.
    
    list flist : list to take indices from
    str/list strings : string to look for, if list looking for several strings
    int nmin, nmax : minimum/maximum index to consider
    bool first : take index of first line meetin condition (not last)
    bool either : if True index lines with any string, otherwise all required
    """
    if nmax is None:
        nmax = len(flist)
    if isinstance(strings, str):
        strings = [strings]
    if either:
        check = any
    else:
        check = all
    for i, item in enumerate(flist[nmin:nmax]):
        if check(string in item for string in strings):
            index = i + nmin
            if first:
                break
    try:
        return index
    except UnboundLocalError:
        UnboundLocalError("Could not find what was sought.")

def strindices(flist, strings, nmin=0, nmax=None, either=False):
    """ Extract indices of all lines in list where strings occur.
    
    list flist : list to take indices from
    str/list strings : string to look for, if list, all must occur in same line
    int nmin, nmax : minimum/maximum index to consider
    bool either : if True index lines with any string, otherwise all required
    """
    if nmax is None:
        nmax = len(flist)
    if isinstance(strings, str):
        strings = [strings]
    if either:
        check = any
    else:
        check = all
    indices = [i+nmin for i, item in enumerate(flist[nmin:nmax]) if
               check(string in item for string in strings)]
    return indices


import ase.io
import sys

class NullDevice():
    """ Blank output stream to redirect the stdout/stderr to """
    def write(self, s):
        pass

def casread(casfile):
    """ ase.io.read() except with no ouput if CASTEP not linked properly """
    oldtargetout = sys.stdout
    oldtargeterr = sys.stderr
    try:
        sys.stdout = NullDevice()
        sys.stderr = NullDevice()
        atoms = ase.io.read(casfile)
        sys.stdout = oldtargetout
        sys.stderr = oldtargeterr
        return atoms
    except:
        """ Bare exception not a problem since error gets raised anyway """
        sys.stdout = oldtargetout
        sys.stderr = oldtargeterr
        raise sys.exc_info()
