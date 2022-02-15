import numpy as np
from misctools import *
from ase import Atoms

recognised_tasks = ['single', 'geometry', "Electronic"]

"""
Module to manage reading of CASTEP input and output files.
Some of the functionality in this module replicates that in ASE, however,
unlike ASE these scripts work with solid solution calculations employing the
virtual crystal approximation (VCA).
"""

############################################################################
#READ CELL FILE
class readcell():
    """ Class for reading CASTEP .cell input files (may have mixed atoms) """
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

    def get_all_calc_params(self):
        """
        Extracts all lines from the .cell file that are NOT related to
        the structure (i.e. not the cell or positions block).
        This is much easier to do for a .cell file (where everything is input)
        than for a .castep file.
        """
        strucblocks = ['%BLOCK LATTICE_', '%BLOCK lattice_',
                       '%block lattice_', '%BLOCK POSITIONS_',
                       '%BLOCK positions_', '%block positions_']
        self.paramlines = []
        ln = 0  # Line number (to iterate through)
        while ln < self.Nlines:
            cellline = self.celllines[ln]
            if all([string not in cellline for string in strucblocks]):
                self.paramlines += [cellline]
                ln += 1
            else:
                ln += 1
                while ('%ENDBLOCK' not in cellline and
                       '%endblock' not in cellline):
                    cellline = self.celllines[ln]
                    ln += 1

    def extract_struc(self, iteration=None):
        """ Ensures behaviour is the same as readcastep """
        return self.casatoms

    def get_kpoints(self):
        """ returns
        list of ints kgrid : k-points per unit cell (MP grid)
        list of floats offset : offset of MP grid
        """
        try:
            lkpts = strindex(self.celllines,
                                  ['kpoints_mp_grid', 'KPOINTS_MP_GRID'],
                                  either=True)
            kgrid = [int(c) for c in self.celllines[lkpts].split()[-3:]]
        except UnboundLocalError:
            kgrid = [8,8,4]  # Defaults to sensible value for RP systems
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
            lpspsb = strindex(self.celllines,
                                   ['%block species_pot',
                                    '%BLOCK species_pot',
                                    '%BLOCK SPECIES_POT'], either=True)
            lpspse = strindex(self.celllines,
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
        lposns = strindex(self.celllines,
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
            mixkey[posstring(self.posns[i, :])] = (elem, {elem: 1.0})
        lposns = strindex(self.celllines,
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
                poskey = posstring(self.posns[i, :])
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
            lpress = strindex(self.celllines,
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
            lconstrs = strindex(self.celllines,
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
#READ CASTEP FILE
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
            lkpts = strindex(self.caslines,
                                  'MP grid size for SCF calculation is')
            kgrid = [int(c) for c in self.caslines[lkpts].split()[-3:]]
        except UnboundLocalError:
            kgrid = [8, 8, 4]
        offset = []
        try:
            loffset = strindex(self.caslines, 'with an offset of')
            offset = [float(o) for o in self.caslines[loffset].split()[-3:]]
        except UnboundLocalError:
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
            lpsps = strindex(self.caslines,
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
        lNions = strindex(self.caslines, 'Total number of ions in cell')
        Nions = int(self.caslines[lNions].split()[7])
        return Nions

    def get_task(self):
        """ returns
        string task : name of task (hopefully one of the recognised_tasks) """
        ltask = strindex(
            self.caslines, 'type of calculation                            :')
        task = self.caslines[ltask].split()[4]
        return task

    def check_complete(self):
        """ returns
        bool complete : True if calculation is complete """
        try:
            strindex(self.caslines, 'Total time          =')
            complete = True
        except UnboundLocalError:
            complete = False
        return complete

    def get_Niterations(self):
        """ returns
        int Niterations : number of structures with enthalpy computed """
        if (self.task == 'single' or self.task == "Electronic"):
            if self.complete == 1:
                Niterations = 1
            else:
                Niterations = 0
        elif self.task == 'geometry':
            lenthalpies = strindices(self.caslines, 'with enthalpy=')
            Niterations = len(lenthalpies)
        return Niterations

    def get_elements(self):
        """ returns
        list of strings : element of each ion (atoms.get_chemical_symbols) """
        lelem = strindex(self.caslines, 'Element ', first=True)
        elems = []
        for casline in self.caslines[lelem+3:lelem+3+self.Nions]:
            elems += [casline.split()[1]]
        return elems

    def get_init_spin(self):
        """ returns
        list of floats spins : initial spin for each ion (in Bohr magnetons)"""
        spins = [0.0]*self.Nions
        try:
            lspin = strindex(self.caslines, 'Initial magnetic')
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
            lspin = strindex(self.caslines,
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

    def get_elempos(self):
        """
        Returns dictionary where key is the element and value if the fractional
        position
        """
        posns = self.get_posns()
        elems = self.get_elements()

        elempos = {}
        for i, pos in enumerate(posns):
            elempos[tuple(pos)] = elems[i]

        return elempos

    def get_mixkey(self, iteration=-1, pureelems=None):
        """ Extract a dictionary mapping mixed atoms onto single site

        int iteration : atom positions (site labels) change during simulation

        purelems: dictionary
           {dopant: pure element} 

        returns
        dict mixkey : mapping -- see mixmap module for more info """
        posns = self.get_posns(iteration=iteration)
        if (self.task == 'single' or self.task == "Electronic"):
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        mixkey = {}
        for i in range(self.Nions):
            elem = self.elems[i]
            # This is the default mixkey for no mixed atoms
            mixkey[posstring(posns[i, :])] = (elem, {elem: 1.0})
        try:
            lmix = strindex(self.caslines, 'Mixture',
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
                            poskey = posstring(posns[i, :])
                    if matchindex is None:
                        raise KeyError('The site ' + posstring(posn) + ' could'
                                       + ' not be matched to any position.')
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
        except (UnboundLocalError, IndexError, TypeError):
            pass  # Normal behaviour if no VCA used

        #Making sure correct pure structure is given
        if pureelems:
            for el in list(pureelems.keys()):
                for pos in posns:
                    key = posstring(pos)
                    if mixkey[key][0] == el:
                        sitemix = mixkey[key][1]
                        mixkey[key] = [pureelems[el], sitemix]

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
                lmax = strindex(self.caslines, 'finished iteration',
                                     first=True, nmin=nmin, nmax=nmax)
                # Above line is to ensure that the geom convergence info
                # is included (occurs after the finished iteration statement)
            else:
                indices = strindices(self.caslines, 'finished iteration',
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
        if (self.task == 'single' or self.task == "Electronic"):
            nmin, nmax = (0, self.Nlines)
        elif self.task == 'geometry':
            nmin, nmax = self.geomrange(iteration=iteration)
        lposn = strindex(self.caslines, 'Element ', nmin=nmin, nmax=nmax,
                first=True)
        poslines = self.caslines[lposn + 3:lposn + 3 + self.Nions]
        for i in range(self.Nions):
            posns[i, :] = [float(p) for p in poslines[i].split()[3:6]]
        return posns

    def get_ext_press(self):
        """ returns
        list of floats press : external pressure in Voigt notation """
        try:
            lpress = strindex(self.caslines,
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
            lconstrs = strindex(self.caslines, 'Cell constraints are:')
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
        if (self.task == 'single' or self.task == "Electronic"):
            nmin, nmax = (0, self.Nlines)
        elif self.task == 'geometry':
            cellconstrs = self.get_cell_constrs()
            if cellconstrs == [0, 0, 0, 0, 0, 0]:
                # at the end or for each iteration since they do not change
                iteration = 0
            nmin, nmax = self.geomrange(iteration=iteration)
        lcell = strindex(self.caslines, 'Real Lattice(A)',
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
            lenthalpy = strindex(self.caslines, 'with enthalpy=',
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
            lenergy = strindex(self.caslines, 'Final energy, E',
                                    nmin=nmin, nmax=nmax)
            energy = float(self.caslines[lenergy].split()[4])
        except UnboundLocalError:
            try:
                lenergy = strindex(self.caslines, 'Final energy =',
                                    nmin=nmin, nmax=nmax)
                energy = float(self.caslines[lenergy].split()[3])
            except UnboundLocalError:
                lenergy = strindex(self.caslines, 'LBFGS: Final Enthalpy')
                energy = float(self.caslines[lenergy].split()[4])

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
        lforce = strindex(self.caslines,
                               ['* Forces *',
                                '* Symmetrised Forces *',
                                '* Constrained Symmetrised Forces *'],
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
        lstress = strindex(self.caslines,
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

    def get_bond_lengths(self):

        """
        This function extracts the bond lengths of the structure after a
        SinglePoint, GeometryOptimization, or Electronic calculation was run.
        Note that CASTEP prints bond lengths only for the optimised structure
        in a geometry optimisation.

        Returns:
        --------
        bondLengths: dict[(atom1, atom2)] = bond length (Å)
            Dictionary whose keys are a tuple of strings containing the
            chemical element and the number in the structure, and the value is
            the bond length, given as a float in Å.
        """
        #Check if calculation completed
        assert self.check_complete(), "The calculation did not complete and" +\
                " the bond lengths cannot be extracted."

        #Get line where bond lengths start being printed
        (_, nmax) = self.geomrange(iteration=-1)
        try:
            #stringstart = "Bond                   Population        Spin       " +\
            #    "Length (A)"
            stringstart='Bond'
            start = strindex(self.caslines, stringstart, nmin=nmax) + 2
        except UnboundLocalError:
            stringstart='Bond'
            #stringstart = "Bond                   Population      Length (A)"
            start = strindex(self.caslines, stringstart, nmin=nmax) + 2

        stringend="Initialisation time ="
        stringend="========================================================="+\
                "============="
        end = strindex(self.caslines, stringend, nmin=nmax+20)

        #Initialise dictionary
        bondLengths = {}

        for line in self.caslines[start:end]:
            cont = line.split()
            key=(cont[0] + cont[1], cont[3] + cont[4])
            bondLengths[key] = float(cont[-1])

        return bondLengths

