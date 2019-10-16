##############################################################################
#### All functions relating to finding the indices at which strings occur ####
##############################################################################

# Extract index of last line in list at which string occurs. Within optional min-max range.
def strindex(list,string,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    i = 0
    for item in list:
        if string in item:
            if i >= nmin and i <= nmax:
                index = i
        i = i+1
    return index

# Extract index of last line in list at which multiple strings occur (all in same line). Within optional min-max range.
def strsindex(list,strings,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    i = 0
    for item in list:
        if all(string in item for string in strings):
            if i >= nmin and i <= nmax:
                index = i
        i = i+1
    return index

# Extract index of last line in list at which multiple strings occur (all in same line). Within optional min-max range.
def strsindex(list,strings,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    i = 0
    for item in list:
        if all(string in item for string in strings):
            if i >= nmin and i <= nmax:
                index = i
        i = i+1
    return index

# Extract index of last line in list at which multiple strings occur (all in same line). Within optional min-max range.
def strsindices(list,strings,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    indices = []
    i = 0
    for item in list:
        if all(string in item for string in strings):
            if i >= nmin and i <= nmax:
                indices.append(i)
        i = i+1
    return indices

# Extract index of last line in list at which one of multiple strings occur. Within optional min-max range.
def eitheror_strsindex(list,strings,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    i = 0
    for item in list:
        if any(string in item for string in strings):
            if i >= nmin and i <= nmax:
                index = i
        i = i+1
    return index

# Extract multiple lists, each with indices of different strings (note, lists possible, if no stypes, assumes alllists!!)
def multilists_strindices(list,strings,nmin=0,nmax=0,shifts=0,stypes=0):
    # Set up lists
    if nmax == 0:
        nmax = len(list)
    Nstrings=len(strings)
    indices = [[]]*Nstrings
    if shifts == 0:
        shifts=[0]*Nstrings
    else:
        if len(shifts) != Nstrings:
            raise ValueError("multilists_strindices: shifts vector not equal length to strings vector.")
    # Each string information: 0=string, 1=anylist, 2=alllist
    if stypes == 0:
        stypes=[(type[s]==list) for s in strings]   # i.e. all lists in strings are alllists
    # Iterate
    i = nmin
    for item in list[nmin:nmax]:
        for s in range(0,Nstrings):
            string=strings[s]
            if stypes[s] == 0:
                if string in item:
                    indices[s]=indices[s]+[i+shifts[s]]              # assuming "string" is a string
            elif stypes[s] == 1:
                if all(substring in item for substring in string):   # assuming "string" is actually a list of strings, want all from list
                    indices[s]=indices[s]+[i+shifts[s]]   
            elif stypes[s] == 2:
                if any(substring in item for substring in string):   # assuming "string" is actually a list of strings, want any from list
                    indices[s]=indices[s]+[i+shifts[s]]   
        i = i+1
    return indices    

# Extract list of indices of lines in list at which string occurs. Within optional min-max range.
def strindices(list,string,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    indices = []
    i = 0
    for item in list:
        if string in item:
            if i >= nmin and i <= nmax:
                indices.append(i)
        i = i+1
    return indices
