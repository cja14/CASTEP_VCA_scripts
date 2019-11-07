"""
Module containing functions to detect the list indices in which strings occur.
Generally the list is normally to be a list of file lines from f.readlines().
"""


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
    return index


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
