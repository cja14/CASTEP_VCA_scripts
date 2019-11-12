import numpy as np
import strindices as stri
from casase import casread
from ase import Atoms
import mixmap

recognised_tasks = ['single', 'geometry']

"""
Module to manage reading of CASTEP input and output files.
Some of the functionality in this module replicates that in ase, however,
unlike ase these scripts work with solid solution calculations employing the
virtual crystal approximation (VCA).
"""


class casatoms(object):
    """ NOTE: treatment of spins has not been implemented in this class. """
    
    def __init__(self, mixatoms=None, mixkey=None, pureatoms=None, spins=None):
        """
        ase.Atoms mixatoms : atomic structure with multiple atoms per site
        dict mixkey : maps between mixed & pure atoms (site or elem format)
        ase.Atoms pureatoms : atomic structure with one atom per site
        list spins : ase.Atoms objects have magmoms attribute (redundant?)
        """
        if spins is not None:
            raise NotImplementedError('Class has not yet been set up to ' +
                                      'work with spins.')
        self.spins = None
        
        if mixkey is not None:
            mkey0 = list(mixkey.keys())[0]
            if len(mkey0.split()) == 3:  # mixkey in site format
                mixkey_type = 'site'
                self.site_mixkey = mixkey
            elif len(mkey0.split()) == 1:  # mixkey in element format
                mixkey_type = 'elem'
            else:
                raise ValueError('Format of mixkey not recognised')
        
        if mixatoms is not None:
            if pureatoms is not None:
                raise KeyError('Cannot specify both mixatoms and pureatoms.')
            self.mixatoms = mixatoms
            if mixkey is None:
                raise KeyError('Cannot specify mixatoms but not mixkey.')
            elif mixkey_type == 'elem':
                self.site_mixkey = mixmap.mixmap.elem2site_mixkey(
                    mixatoms, mixkey, pure=False)
        
        elif pureatoms is not None:
            self.pureatoms = pureatoms
            if mixkey is None:
                self.site_mixkey = mixmap.mixmap.trivial_mixkey(pureatoms)
                self.mixatoms = pureatoms
            else:
                self.mixatoms = mixmap.create_mixture(pureatoms, mixkey)
            if mixkey_type == 'elem':
                self.site_mixkey = mixmap.mixmap.elem2site_mixkey(
                    pureatoms, mixkey, pure=True)
        
        else:
            raise KeyError('Both mixatoms and pureatoms cannot be None.')
        
        self.mapping = mixmap.mixmap(self.mixatoms, self.site_mixkey)
        try:
            self.elem_mixkey = mixmap.mixmap.site2elem_mixkey(self.site_mixkey)
        except ValueError:
            self.elem_mixkey = None
        self.pureatoms = self.mapping.mix2pure(self.mixatoms)
        self.paramlines = []
        self.celllines = None
        self.Nlines = None
    
    def __repr__(self):
        """ return pureatoms and mixkey? """
        if self.elem_mixkey is not None:
            mixkey = self.elem_mixkey
        else:
            mixkey = self.site_mixkey
        return self.__class__.__name__ + '(' + self.pureatoms.__repr__() + \
            ', mixkey=' + str(mixkey) + ')'
    
    def update_mix(self):
        """ Update mixatoms when changes have been made to pureatoms. """
        self.mixatoms = self.mapping.pure2mix(self.pureatoms)
    
    def update_pure(self):
        """ Update pureatoms when changes have been made to mixatoms. """
        self.pureatoms = self.mapping.mix2pure(self.mixatoms)
    
    def get_mixatoms(self):
        return self.mixatoms
    
    def set_mixatoms(self, mixatoms):
        self.mixatoms = mixatoms
        self.update_pure()
    
    def get_pureatoms(self):
        return self.pureatoms
    
    def set_pureatoms(self, pureatoms):
        self.pureatoms = pureatoms
        self.update_mix()
        
    def extract_cellparams(self):
        """ Extracts all lines from the .cell file that are NOT related to
        the structure (i.e. not the cell or positions block).
        
        This is much easier to do for a .cell file (where everything is input)
        than for a .castep file."""
        strucblocks = ['%block lattice_', '%block positions_']
        self.paramlines = []
        ln = 0  # Line number (to iterate through)
        while ln < self.Nlines:
            cellline = self.celllines[ln]
            ln += 1
            if all([string not in cellline.lower() for string in strucblocks]):
                self.paramlines += [cellline]
            else:
                while '%endblock' not in cellline.lower():
                    cellline = self.celllines[ln]
                    ln += 1
    
    def read_cellparams(self, cellfile):
        """ Wraps the extract_cellparams method when a .cell file not open """
        self.celllines = open(cellfile, 'r').readlines()
        self.Nlines = len(self.celllines)
        self.extract_cellparams()
        self.celllines = None
        self.Nlines = None
        
    def get_keyword(self, keyword, value=False):
        """ Extract from paramlines the value associated with a keyword.
        If keyword is present (but has no value) method will return True.
        If keyword is not present, returns optional input value"""
        try:
            idx = stri.strindex(self.paramlines, keyword)
            ksplit = [s for s in self.paramlines[idx].split() if s not in
                      ['=', ':']]
            if len(ksplit) == 1:
                value = True
            elif len(ksplit) == 2:
                value = ksplit[1]
            else:
                value = ksplit[1:]
        except (IndexError, UnboundLocalError):
            pass
        return value
    
    def set_keyword(self, keyword, value=True):
        """ Set the value associated with a keyword.
        If value is not specified as optional input, the keyword will be
        present in paramlines but with no associated value."""
        try:
            idx = stri.strindex(self.paramlines, keyword)
        except (IndexError, UnboundLocalError):
            idx = len(self.paramlines)
            self.paramlines += ['\n'*bool(value), '\n'*bool(value)]
        if value:
            if isinstance(value, bool):
                self.paramlines[idx] = keyword + '\n'
            else:
                self.paramlines[idx] = keyword + ' = ' + str(value) + '\n'
        else:
            end = idx + 1
            if self.paramlines[end:] and not self.paramlines[end].split():
                end = idx + 2
            self.paramlines = self.paramlines[:idx] + self.paramlines[end:]
    
    def get_block(self, keyword, value=None):
        try:
            idxs = stri.strindices(self.paramlines, keyword)
            start, end = idxs[0], idxs[1]
            value = []
            for idx in range(start+1, end):
                value += [self.paramlines[idx].split()]
        except (IndexError, UnboundLocalError):
            pass
        return value
    
    def set_block(self, keyword, values):
        try:
            idxs = stri.strindices(self.paramlines, keyword)
            start, end = idxs[0], idxs[1]
            before = self.paramlines[:start]
            after = self.paramlines[end+1:]
        except (IndexError, UnboundLocalError):
            before = self.paramlines[:]
            after = ['\n'*bool(values)]
        if values:
            block = ['%BLOCK ' + keyword + '\n']
            for val in values:
                if isinstance(val, str):
                    block += [val + '\n']
                elif isinstance(val, int) or isinstance(val, float):
                    block += [str(val) + '\n']
                else:
                    block += [' '.join([str(v) for v in val]) + '\n']
            block += ['%ENDBLOCK ' + keyword + '\n']
        else:
            block = []
            if after and not after[0].split():
                after = after[1:]
        self.paramlines = before + block + after
    
    def writecell(self, cellfile, pure=False, frac=True):
        """ Write a CASTEP .cell input file """
        # Print cell
        caslines = ['%BLOCK LATTICE_CART\n']
        cell = self.mixatoms.get_cell()
        for i in range(3):
            caslines += ['\t'+'\t'.join([str('{0:.8f}'.format(cell[i, j]))
                                         for j in range(3)])+'\n']
        caslines += ['%ENDBLOCK LATTICE_CART\n']+['\n']
        if pure:
            atoms = self.pureatoms
        else:
            atoms = self.mixatoms
        elems = atoms.get_chemical_symbols()
        if frac:
            caslines += ['%BLOCK POSITIONS_FRAC\n']
            mixposns = atoms.get_scaled_positions()
        else:
            caslines += ['%BLOCK POSITIONS_ABS\n']
            mixposns = atoms.get_positions()
        for i, elem in enumerate(elems):
            (m, wt) = self.mapping.mixsitemixes[i]
            if pure:
                wt = 1.0
            if wt == 1.0:
                mixstring = ''
            elif wt <= 0.0 or wt > 1.0:
                raise ValueError('Trying to print ion index ' + str(i) +
                                 ' with weight ' + str(wt) +
                                 ' and do not know what to do.')
            else:
                mixstring = '\tMIXTURE=('+str(m)+' '+str(wt)+')'
            if self.spins is not None:
                spin = self.spins[i]
            else:
                spin = 0
            if spin == 0:
                spinstring = ''
            else:
                spinstring = '\tSPIN='+str(spin)
            posstring = '\t'.join(
                [str('{0:.8f}'.format(mixposns[i, j])) for j in range(3)])
            caslines += ['\t' + elem + '\t' + posstring + spinstring +
                         mixstring + '\n']
        if frac:
            caslines += ['%ENDBLOCK POSITIONS_FRAC\n']
        else:
            caslines += ['%ENDBLOCK POSITIONS_ABS\n']
        caslines += self.paramlines
        open(cellfile, 'w').writelines(caslines)


class readcell(casatoms):
    """ Class for reading CASTEP .cell files (may have mixed atoms) """
    def __init__(self, cellfile, flttol=1e-4):
        """
        string cellfile : path to .cell file
        float flttol : two numbers considered equal within this tolerance
        """
        self.celllines = open(cellfile, 'r').readlines()
        self.Nlines = len(self.celllines)
        self.flttol = flttol
        self.casatoms = casread(cellfile)
        self.elems = self.casatoms.get_chemical_symbols()
        self.posns = self.casatoms.get_scaled_positions()
        self.Nions = len(self.casatoms)
        self.get_all_calc_params()
        
    def extract_struc(self, iteration=None):
        """ Ensures behaviour is the same as readcastep """
        return self.casatoms
    
    def get_kpoints(self):
        """ returns
        list of ints kgrid : k-points per unit cell (MP grid)
        list of floats offset : offset of MP grid
        """
        try:
            lkpts = stri.strindex(self.celllines,
                                  ['kpoints_mp_grid', 'KPOINTS_MP_GRID'],
                                  either=True)
            kgrid = [int(c) for c in self.celllines[lkpts].split()[-3:]]
        except UnboundLocalError:
            kgrid = [5, 5, 1]  # Defaults to sensible value for RP systems
        offset = []
        for k in kgrid:
            if k % 2 == 0:
                offset += [0.0]
            else:
                offset += [1.0/(2*k)]
        return kgrid, offset
    
    def get_psps(self):
        """ returns
        list of strings pseudos : CASTEP pseudo-potential strings """
        try:
            lpspsb = stri.strindex(self.celllines,
                                   ['%block species_pot',
                                    '%BLOCK species_pot',
                                    '%BLOCK SPECIES_POT'], either=True)
            lpspse = stri.strindex(self.celllines,
                                   ['%endblock species_pot',
                                    '%ENDBLOCK species_pot',
                                    '%ENDBLOCK SPECIES_POT'], either=True)
            pseudos = {}
            for i in range(lpspsb+1, lpspse):
                elem, psp = tuple(self.celllines[i].split())
                pseudos[elem] = psp
        except UnboundLocalError:
            pseudos = None
        return pseudos
    
    def get_elements(self):
        """ returns
        list of strings : element of each ion (atoms.get_chemical_symbols) """
        return self.elems
    
    def get_init_spin(self):
        """ returns
        list of floats spins : initial spin for each ion (in Bohr magnetons)"""
        spins = [0.0]*self.Nions
        lposns = stri.strindex(self.celllines,
                               ['%block positions_frac',
                                '%BLOCK positions_frac',
                                '%BLOCK POSITIONS_FRAC',
                                '%block positions_abs',
                                '%BLOCK positions_abs',
                                '%BLOCK POSITIONS_ABS'], either=True)
        for i in range(self.Nions):
            if (len(self.celllines[lposns+1+i].split()) > 4 and
                ('SPIN' in self.celllines[lposns+1+i] or
                 'spin' in self.celllines[lposns+1+i])):
                lnsplt = self.celllines[lposns+1+i].split()
                for j, string in enumerate(lnsplt):
                    if 'SPIN' in string or 'spin' in string:
                        break
                spins[i] = float(lnsplt[j].split('=')[1])
        return spins
    
    def get_mixkey(self, iteration=None):
        """ Extract a dictionary mapping mixed atoms onto single site
        returns
        dict mixkey : mapping -- see mixmap module for more info """
        mixkey = {}
        for i in range(self.Nions):
            elem = self.elems[i]
            # This is the default mixkey for no mixed atoms
            mixkey[mixmap.posstring(self.posns[i, :])] = (elem, {elem: 1.0})
        lposns = stri.strindex(self.celllines,
                               ['%block positions_frac',
                                '%BLOCK positions_frac',
                                '%BLOCK POSITIONS_FRAC',
                                '%block positions_abs',
                                '%BLOCK positions_abs',
                                '%BLOCK POSITIONS_ABS'], either=True)
        for i in range(self.Nions):
            if (len(self.celllines[lposns+1+i].split()) > 4 and
                ('MIXTURE' in self.celllines[lposns+1+i] or
                 'mixture' in self.celllines[lposns+1+i])):
                lnsplt = self.celllines[lposns+1+i].split()
                for j, string in enumerate(lnsplt):
                    if 'MIXTURE' in string or 'mixture' in string:
                        imix = j
                elem = self.elems[i]
                wt = float(lnsplt[imix+1].replace(')', ''))
                poskey = mixmap.posstring(self.posns[i, :])
                siteelem, wts = mixkey[poskey]
                # siteelem could be overwritten when sorting
                wts[elem] = wt
                elemkey = sorted(list(set(wts.keys())))[0]
                mixkey[poskey] = (elemkey, wts)
        return mixkey
    
    def get_posns(self, iteration=-1):
        """ Takes iteration so compatable with readcastep
        returns
        np.array(Nions, 3) posns : fractional position of each ion """
        return self.posns
    
    def get_ext_press(self):
        """ returns
        list of floats press : external pressure in Voigt notation """
        try:
            lpress = stri.strindex(self.celllines,
                                   ['%block external_pressure',
                                    '%BLOCK external_pressure',
                                    '%BLOCK EXTERNAL_PRESSURE'], either=True)
            if len(self.celllines[lpress+1].split()) == 1:
                lpress += 1
            presslines = self.celllines[lpress+1:lpress+4]
            press = [0.0]*6
            press[0] = float(presslines[0].split()[0])
            press[1] = float(presslines[1].split()[0])
            press[2] = float(presslines[2].split()[0])
            press[3] = float(presslines[1].split()[1])
            press[4] = float(presslines[0].split()[2])
            press[5] = float(presslines[0].split()[1])
        except UnboundLocalError:
            press = [0.0]*6
        return press
    
    def get_cell_constrs(self):
        """ returns
        list of ints cellconstrs : CASTEP cell constraints (0 = fixed) """
        try:
            lconstrs = stri.strindex(self.celllines,
                                     ['%block cell_constraints',
                                      '%BLOCK cell_constraints',
                                      '%BLOCK CELL_CONSTRAINTS'], either=True)
            cellconstrs = [int(c) for c in self.celllines[lconstrs+1].split()]
            cellconstrs += [int(c) for c in self.celllines[lconstrs+2].split()]
        except UnboundLocalError:
            cellconstrs = [1, 2, 3, 4, 5, 6]  # Equates to no constraints
        return cellconstrs

    def get_cell(self, iteration=-1):
        """ returns
        np.array(3, 3) cell : unit cell vectors (Angstroms) """
        return self.casatoms.get_cell()

#########################################################################


class readcas():
    """
    Class for extracting info from .castep output files (may have mixed atoms)
    """
    
    def __init__(self, casfile, flttol=1e-4):
        """
        string cellfile : path to .castep file
        float flttol : two numbers considered equal within this tolerance
        """
        self.caslines = open(casfile, 'r').readlines()
        self.Nlines = len(self.caslines)
        self.Nions = self.get_Nions()
        self.task = self.get_task()
        self.flttol = flttol  # For comparing floats
        if self.task not in recognised_tasks:
            raise ValueError('Do not recognise task:' + self.task)
        self.complete = self.check_complete()
        self.Niterations = self.get_Niterations()
        self.elems = self.get_elements()
    
    def extract_struc(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        ase.Atoms casatoms : structure at the desired iteration """
        posns = self.get_posns(iteration=iteration)
        cell = self.get_cell(iteration=iteration)
        casatoms = Atoms(scaled_positions=posns, cell=cell,
                         symbols=self.elems, pbc=True)
        return casatoms
    
    def get_kpoints(self):
        """ returns
        list of ints kgrid : k-points per unit cell (MP grid)
        list of floats offset : offset of MP grid """
        try:
            lkpts = stri.strindex(self.caslines,
                                  'MP grid size for SCF calculation is')
            kgrid = [int(c) for c in self.caslines[lkpts].split()[-3:]]
        except UnboundLocalError:
            kgrid = [5, 5, 1]
        offset = []
        for k in kgrid:
            if k % 2 == 0:
                offset += [0.0]
            else:
                offset += [1.0/(2*k)]
        return kgrid, offset

    def get_psps(self):
        """ returns
        list of strings pseudos : CASTEP pseudo-potential strings """
        try:
            lpsps = stri.strindex(self.caslines,
                                  'Files used for pseudopotentials:')
            pseudos = {}
            i = lpsps+1
            while self.caslines[i].split():
                elem, psp = tuple(self.caslines[i].split())
                pseudos[elem] = psp
                i += 1
        except UnboundLocalError:
            pseudos = None
        return pseudos
    
    def get_Nions(self):
        """ returns
        int Nions : number of ions in cell """
        lNions = stri.strindex(self.caslines, 'Total number of ions in cell')
        Nions = int(self.caslines[lNions].split()[7])
        return Nions

    def get_task(self):
        """ returns
        string task : name of task (hopefully one of the recognised_tasks) """
        ltask = stri.strindex(
            self.caslines, 'type of calculation                            :')
        task = self.caslines[ltask].split()[4]
        return task

    def check_complete(self):
        """ returns
        bool complete : True if calculation is complete """
        try:
            stri.strindex(self.caslines, 'Total time          =')
            complete = True
        except UnboundLocalError:
            complete = False
        return complete
    
    def get_Niterations(self):
        """ returns
        int Niterations : number of structures with enthalpy computed """
        if self.task == 'single':
            if self.complete == 1:
                Niterations = 1
            else:
                Niterations = 0
        elif self.task == 'geometry':
            lenthalpies = stri.strindices(self.caslines, 'with enthalpy=')
            Niterations = len(lenthalpies)
        return Niterations
    
    def get_elements(self):
        """ returns
        list of strings : element of each ion (atoms.get_chemical_symbols) """
        lelem = stri.strindex(self.caslines, 'Element ', first=True)
        elems = []
        for casline in self.caslines[lelem+3:lelem+3+self.Nions]:
            elems += [casline.split()[1]]
        return elems
    
    def get_init_spin(self):
        """ returns
        list of floats spins : initial spin for each ion (in Bohr magnetons)"""
        spins = [0.0]*self.Nions
        try:
            lspin = stri.strindex(self.caslines, 'Initial magnetic')
            spinlines = self.caslines[lspin+3:lspin+3+self.Nions]
            for i in range(self.Nions):
                spins[i] = float(spinlines[i].split()[4])
        except UnboundLocalError:
            pass  # if no initial spins this table won't appear
        return spins
    
    def get_final_spin(self):
        """ returns
        list of floats spins : final spin for each ion (in Bohr magnetons) """
        spins = [0.0]*self.Nions
        try:
            lspin = stri.strindex(self.caslines,
                                  'Atomic Populations (Mulliken)')
            if (('spin' in self.caslines[lspin+2] or 'Spin'
                 in self.caslines[lspin+2])):
                spinlines = self.caslines[lspin+4:lspin+4+self.Nions]
                strfactor = self.caslines[lspin+2].split()[-1]
                if strfactor == '(hbar)':
                    fltfactor = 2.0
                elif strfactor == '(hbar/2)':
                    fltfactor = 1.0
                else:
                    raise ValueError('Scale factor for spins: ' + strfactor +
                                     ' not recognised.')
                for i in range(self.Nions):
                    spins[i] = float(spinlines[i].split()[-1])*fltfactor
        except UnboundLocalError:
            raise UnboundLocalError(
                'Could not find final atomic populations,' +
                ' are you sure the calcation completed?')
        return spins
    
    def get_mixkey(self, iteration=-1):
        """ Extract a dictionary mapping mixed atoms onto single site
        
        int iteration : atom positions (site labels) change during simulation
        
        returns
        dict mixkey : mapping -- see mixmap module for more info """
        posns = self.get_posns(iteration=iteration)
        if self.task == 'single':
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        mixkey = {}
        for i in range(self.Nions):
            elem = self.elems[i]
            # This is the default mixkey for no mixed atoms
            mixkey[mixmap.posstring(posns[i, :])] = (elem, {elem: 1.0})
        try:
            lmix = stri.strindex(self.caslines, 'Mixture',
                                 nmin=nmin, nmax=nmax)
            l = lmix + 3  # l is line index (which we'll iterate through)
            wts = {}
            matchindex = None
            posn = None
            # Whilst in mixture block
            while self.caslines[l].split()[0] == 'x':
                mixline = self.caslines[l].split()
                if len(mixline) == 8:
                    if posn is not None:  # Then poskey must be defined
                        mixkey[poskey] = (self.elems[matchindex], wts)
                        posn, matchindex = None, None
                        wts = {}
                    posn = [float(p) for p in mixline[2:5]]
                    elem = mixline[5]
                    for i in range(self.Nions):
                        dist = np.linalg.norm(np.array(posn) - posns[i, :])
                        if (self.elems[i] == elem and dist < self.flttol):
                            matchindex = i
                            # Since it is the position from posns the key
                            # would be written for
                            poskey = mixmap.posstring(posns[i, :])
                    if matchindex is None:
                        raise KeyError('The site ' + mixmap.posstring(posn) +
                                       " couldn't be matched to any position.")
                    wt = mixline[6]
                else:
                    elem = mixline[1]
                    wt = mixline[2]
                wts[elem] = float(wt)
                l += 1
            # Stopped reading the file but the last mix is probably still open
            if posn is not None:
                mixkey[poskey] = (self.elems[matchindex], wts)
                posn, matchindex = None, None
                wts = {}
        except (UnboundLocalError, IndexError):
            pass  # Normal behaviour if no VCA used
        return mixkey
    
    def geomrange(self, iteration=-1, nmin=0, nmax=None):
        """
        int iteration : positive count from front and negative from back
        Note: iteration == None means an unconstrained data extraction
        (includes hanging/incomplete geom iterations!)
        int nmin, nmax : minimum/maximum line indices to consider
        Note: will count WITHIN these indices!
        
        returns
        int (lmin, lmax) : minimum/maximum line index for iteration """
        if nmax is None:
            nmax = self.Nlines
        if iteration is None:
            lmin = nmin
            lmax = nmax
        else:
            if ((abs(iteration) > self.Niterations or
                 iteration == self.Niterations)):
                raise IndexError('Cannot extract information for iteration '
                                 + str(iteration) + ' since only ' +
                                 str(self.Niterations) +
                                 ' have been performed.')
            if iteration == 0 or iteration == -self.Niterations:
                lmin = nmin
                lmax = stri.strindex(self.caslines, 'finished iteration',
                                     first=True, nmin=nmin, nmax=nmax)
                # Above line is to ensure that the geom convergence info
                # is included (occurs after the finished iteration statement)
            else:
                indices = stri.strindices(self.caslines, 'finished iteration',
                                          nmin=nmin, nmax=nmax)
                lmin = indices[iteration - 1]
                lmax = indices[iteration]
        return (lmin, lmax)
    
    def get_posns(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        returns
        np.array(Nions, 3) posns : fractional position of each ion """
        posns = np.zeros((self.Nions, 3))
        if self.task == 'single':
            nmin, nmax = (0, self.Nlines)
        elif self.task == 'geometry':
            nmin, nmax = self.geomrange(iteration=iteration)
        lposn = stri.strindex(self.caslines, 'Element ', nmin=nmin, nmax=nmax)
        poslines = self.caslines[lposn + 3:lposn + 3 + self.Nions]
        for i in range(self.Nions):
            posns[i, :] = [float(p) for p in poslines[i].split()[3:6]]
        return posns
    
    def get_ext_press(self):
        """ returns
        list of floats press : external pressure in Voigt notation """
        try:
            lpress = stri.strindex(self.caslines,
                                   'External pressure/stress (GPa)')
            presslines = self.caslines[lpress+1:lpress+4]
            press = [0.0]*6
            press[0] = float(presslines[0].split()[0])
            press[1] = float(presslines[1].split()[0])
            press[2] = float(presslines[2].split()[0])
            press[3] = float(presslines[1].split()[1])
            press[4] = float(presslines[0].split()[2])
            press[5] = float(presslines[0].split()[1])
        except UnboundLocalError:
            press = [0.0]*6
        return press
    
    def get_cell_constrs(self):
        """ returns
        list of ints cellconstrs : CASTEP cell constraints (0 = fixed) """
        try:
            lconstrs = stri.strindex(self.caslines, 'Cell constraints are:')
            cellconstrs = [int(c) for c in
                           self.caslines[lconstrs].split()[3:9]]
        except UnboundLocalError:
            cellconstrs = None
        return cellconstrs
    
    def get_cell(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        np.array(3, 3) cell : unit cell vectors (Angstroms) """
        cell = np.zeros((3, 3))
        if self.task == 'single':
            nmin, nmax = (0, self.Nlines)
        elif self.task == 'geometry':
            cellconstrs = self.get_cell_constrs()
            if cellconstrs == [0, 0, 0, 0, 0, 0]:
                # at the end or for each iteration since they do not change
                iteration = 0
            nmin, nmax = self.geomrange(iteration=iteration)
        lcell = stri.strindex(self.caslines, 'Real Lattice(A)',
                              nmin=nmin, nmax=nmax)
        celllines = self.caslines[lcell + 1:lcell + 4]
        for i in range(3):
            cell[i, :] = [float(p) for p in celllines[i].split()[0:3]]
        return cell
    
    def get_enthalpy(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        float enthalpy : cell enthalpy in eV """
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract enthalpy.')
        if self.task == 'single':
            raise NotImplementedError(
                'get_enthalpy not implemented if task == single.')
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
            lenthalpy = stri.strindex(self.caslines, 'with enthalpy=',
                                      nmin=nmin, nmax=nmax)
            enthalpy = float(self.caslines[lenthalpy].split()[6])
        return enthalpy

    def get_energy(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        float energy : cell energy in eV """
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract enthalpy.')
        if self.task == 'single':
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        try:
            lenergy = stri.strindex(self.caslines, 'Final energy, E',
                                    nmin=nmin, nmax=nmax)
            energy = float(self.caslines[lenergy].split()[4])
        except UnboundLocalError:
            lenergy = stri.strindex(self.caslines, 'Final energy =',
                                    nmin=nmin, nmax=nmax)
            energy = float(self.caslines[lenergy].split()[3])
        return energy
    
    def get_forces(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        np.array(Nions, 3) forces : force vector of each ion (eV/Ang) """
        forces = np.zeros((self.Nions, 3))
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract forces.')
        if self.task == 'single':
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        lforce = stri.strindex(self.caslines,
                               ['* Forces *',
                                '* Symmetrised Forces *'],
                               either=True, nmin=nmin, nmax=nmax)
        forcelines = self.caslines[lforce+6:lforce+6+self.Nions]
        for i in range(self.Nions):
            newforcelinesplt = [f for f in forcelines[i].split()
                                if f != '(mixed)']
            forces[i, :] = [float(p.replace('(cons\'d)', ''))
                            for p in newforcelinesplt[3:6]]
        return forces

    def get_stresses(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        np.array(3, 3) stresses : stress matrix (eV/Ang^3) """
        stresses = np.zeros((3, 3))
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract stresses.')
        if self.task == 'single':
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        lstress = stri.strindex(self.caslines,
                                ['* Stress Tensor *',
                                 '* Symmetrised Stress Tensor *'],
                                either=True, nmin=nmin, nmax=nmax)
        stresslines = self.caslines[lstress+6:lstress+9]
        for i in range(3):
            stresses[i, :] = [float(p) for p in stresslines[i].split()[2:5]]
        return stresses
    
    def get_Fmax(self, iteration=-1):
        """
        int iteration : index of desired iteration in simulation
        
        returns
        float Fmax : maximum force on any ion (eV/Ang) """
        forces = self.get_forces(iteration=iteration)
        return max(np.linalg.norm(forces, axis=1))
