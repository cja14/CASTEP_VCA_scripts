import numpy as np
import strindices as stri
from casase import casread

recognised_tasks = ['single', 'geometry']


def pzero(x):
    """ Make sure zeros are displayed as positive """
    if x == 0.0:
        x = 0.0
    return x


def posstring(posn):
    """ unambiguously flattern a position array to a string """
    return ' '.join([str('{0:.6f}'.format(pzero(posn[j])))
                     for j in range(len(posn))])


class readcell():
    def __init__(self, cellfile, flttol=1e-4):
        self.celllines = open(cellfile, 'r').readlines()
        self.Nlines = len(self.celllines)
        self.flttol = flttol
        self.casatoms = casread(cellfile)
        self.elems = self.casatoms.get_chemical_symbols()
        self.posns = self.casatoms.get_scaled_positions()
        self.Nions = len(self.casatoms)

    def extract_struc(self, iteration=None):
        return self.casatoms

    def get_kpoints(self):
        try:
            lkpts = stri.eitheror_strsindex(self.celllines,
                                            ['kpoints_mp_grid',
                                             'KPOINTS_MP_GRID'])
            kgrid = [int(c) for c in self.celllines[lkpts].split()[-3:]]
        except UnboundLocalError:
            kgrid = [5, 5, 1]
        offset = []
        for k in kgrid:
            if k % 2 == 0:
                offset += [0.0]
            else:
                offset += [1.0/(2*k)]
        return kgrid,offset

    def get_psps(self):
        try:
            lpspsb = stri.eitheror_strsindex(self.celllines,
                                             ['%block species_pot',
                                              '%BLOCK species_pot',
                                              '%BLOCK SPECIES_POT'])
            lpspse = stri.eitheror_strsindex(self.celllines,
                                             ['%endblock species_pot',
                                              '%ENDBLOCK species_pot',
                                              '%ENDBLOCK SPECIES_POT'])
            pseudos = {}
            for i in range(lpspsb+1, lpspse):
                elem, psp = tuple(self.celllines[i].split())
                pseudos[elem] = psp
        except UnboundLocalError:
            pseudos = None
        return pseudos
    
    def get_elements(self):
        return self.elems

    # Get initial magnetic moment of each atom
    def get_init_spin(self):
        #raise NotImplementedError('Does not yet work for spins.')
        return [0]*self.Nions
    
    # Extract a mixkey for every ionic site in the unitcell
    def get_mixkey(self, iteration=None):
        mixkey = {}
        for i in range(self.Nions):
            elem = self.elems[i]
            # This is the default mixkey for no mixed atoms
            mixkey[posstring(self.posns[i, :])] = {elem: {elem: 1.0}}
        lposns = stri.eitheror_strsindex(self.celllines,
                                         ['%block positions_frac',
                                          '%BLOCK positions_frac',
                                          '%BLOCK POSITIONS_FRAC',
                                          '%block positions_abs',
                                          '%BLOCK positions_abs',
                                          '%BLOCK POSITIONS_ABS'])
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
                fullmatchdict = mixkey[poskey]
                matchdict = list(fullmatchdict.values())[0]
                matchdict[elem] = wt
                elemkey = sorted(list(set(matchdict.keys())))[0]
                mixkey[poskey] = {elemkey: matchdict}
        return mixkey
    
    def get_posns(self, iteration=-1):
        return self.posns

    # Get the external cell pressure used in the calculation
    def get_ext_press(self):
        try:
            lpress = stri.eitheror_strsindex(self.celllines,
                                           ['%block external_pressure',
                                            '%BLOCK external_pressure',
                                            '%BLOCK EXTERNAL_PRESSURE'])
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
        except:
            press = [0.0]*6
        return press
    
    # Get the Cell Constraints for a Geometry Optimisation
    def get_cell_constrs(self):
        try:
            lconstrs = stri.eitheror_strsindex(self.celllines,
                                               ['%block cell_constraints',
                                                '%BLOCK cell_constraints',
                                                '%BLOCK CELL_CONSTRAINTS'])
            cellconstrs = [int(c) for c in self.celllines[lconstrs+1].split()]
            cellconstrs += [int(c) for c in self.celllines[lconstrs+2].split()]
        except UnboundLocalError:
            cellconstrs = [1, 2, 3, 4, 5, 6]
        return cellconstrs

    def get_cell(self, iteration=-1):
        return self.casatoms.get_cell()

#########################################################################


class readcas():
    """
    Class for extracting info from .castep output files (may have mixed atoms)
    """
    def __init__(self, casfile, flttol=1e-4):
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
        # so can use this class without loading the module import ase here
        from ase import Atoms
        posns = self.get_posns(iteration=iteration)
        cell = self.get_cell(iteration=iteration)
        casatoms = Atoms(scaled_positions=posns, cell=cell,
                         symbols=self.elems, pbc=True)
        return casatoms

    def get_kpoints(self):
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
        lNions = stri.strindex(self.caslines, 'Total number of ions in cell')
        N = int(self.caslines[lNions].split()[7])
        return N

    def get_task(self):
        ltask = stri.strindex(
            self.caslines, 'type of calculation                            :')
        task = self.caslines[ltask].split()[4]
        return task

    def check_complete(self):
        try:
            stri.strindex(self.caslines, 'Total time          =')
            complete = True
        except UnboundLocalError:
            complete = False
        return complete
    
    # Basically just count number of enthalpies generated
    def get_Niterations(self):
        if self.task == 'single':
            if self.complete == 1:
                N = 1
            else:
                N = 0
        elif self.task == 'geometry':
            lenthalpies = stri.strindices(self.caslines, 'with enthalpy=')
            N = len(lenthalpies)
        return N
    
    def get_elements(self):
        lelem = stri.strindices(self.caslines, 'Element')[0]
        elems = []
        for casline in self.caslines[lelem+3:lelem+3+self.Nions]:
            elems += [casline.split()[1]]
        return elems

    # Get initial magnetic moment of each atom
    def get_init_spin(self):
        spins = [0.0]*self.Nions
        try:
            lspin = stri.strindex(self.caslines, 'Initial magnetic')
            spinlines = self.caslines[lspin+3:lspin+3+self.Nions]
            for i in range(self.Nions):
                spins[i] = float(spinlines[i].split()[4])
        except UnboundLocalError:
            pass  # if no initial spins this table won't appear
        return spins

    # Get final magnetic moment (in Bohr magnetons) of each atom
    def get_final_spin(self):
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
    
    # Make sure zeros are displayed as positive
    def pzero(self, x):
        if x == 0.0:
            x = 0.0
        return x
        
    # unambiguously flattern a position array to a string
    def posstring(self, posn):
        return ' '.join([str(pzero(posn[j])) for j in range(len(posn))])
    
    # Extract a mixkey for every ionic site in the unitcell
    # Actually the iteraction is important since the posns of ions will change
    def get_mixkey(self, iteration=-1):
        posns = self.get_posns(iteration=iteration)
        if self.task == 'single':
            (nmin, nmax) = (0, self.Nlines)
        elif self.task == 'geometry':
            (nmin, nmax) = self.geomrange(iteration=iteration)
        mixkey = {}
        for i in range(self.Nions):
            elem = self.elems[i]
            # This is the default mixkey for no mixed atoms
            mixkey[posstring(posns[i, :])] = {elem: {elem: 1.0}}
        try:
            lmix = stri.strindices(self.caslines,
                                   'Mixture', nmin=nmin, nmax=nmax)[-1]
            l = lmix + 3
            wts = {}
            matchindex = None
            posn = None
            while self.caslines[l].split()[0] == 'x':
                mixline = self.caslines[l].split()
                if len(mixline) == 8:
                    if posn is not None:
                        mixkey[poskey] = {self.elems[matchindex]: wts}
                        # This step isn't strictly necessary but to be safe
                        posn, matchindex = None, None
                        wts = {}
                    # Must be a list (cannot have dict of np array type)
                    posn = [float(p) for p in mixline[2:5]]
                    elem = mixline[5]
                    # It seems ridiculous but in the same .castep file
                    # the same number can be printed to multiple
                    # different levels of accuracy therefore a floating
                    # point tolerance is still required
                    for i in range(self.Nions):
                        if (self.elems[i] == elem and
                            all([abs(posn[j]-posns[i, j]) < self.flttol
                                 for j in range(3)])):
                            matchindex = i
                            # Since it is the position from posns the key
                            # would be written for
                            poskey = posstring(posns[i, :])
                    if matchindex is None:
                        raise KeyError('The site ' + posstring(posn) +
                                       ' could not be matched to any position.')
                    wt = mixline[6]
                else:
                    elem = mixline[1]
                    wt = mixline[2]
                wts[elem] = float(wt)
                l += 1
            # Stopped reading the file but the last mix is probably still open
            if posn is not None:
                mixkey[poskey] = {self.elems[matchindex]: wts}
                # This step isn't strictly necessary but just to be safe
                posn,matchindex = None,None
                wts = {}                
        except IndexError:
            pass  # This isn't a problem, in fact it's normal behaviour
        return mixkey
            
    # Find geom iteration range: iterations now work as you would expect
    # => negative indices count back from end!
    # NOTE: iteration=None means an unconstrained data extraction
    # (includes hanging/incomplete geom iterations!)
    def geomrange(self, iteration=-1, nmin=0, nmax=0):
        if nmax == 0:
            nmax = self.Nlines
        # WARNING: unconstrained data extractions can pick up the last
        # instance of that data from anywhere in the file
        # therefore different parameters are not guaranteed to be
        # consistent with one another!
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
                lmax = stri.strindices(self.caslines, 'finished iteration')[0]
                # Above line is to ensure that the geom convergence info
                # is included (occurs after the finished iteration statement)
            else:
                lmin = stri.strindices(self.caslines,
                                       'finished iteration')[iteration - 1]
                lmax = stri.strindices(self.caslines,
                                       'finished iteration')[iteration]
        return (lmin, lmax)

    def get_posns(self,iteration=-1):
        posns = np.zeros((self.Nions, 3))
        if self.task == 'single':
            nmin, nmax = (0, self.Nlines)
        elif self.task == 'geometry':
            nmin, nmax = self.geomrange(iteration=iteration)
        lposn = stri.strindex(
            self.caslines,
            'Element    Atom        Fractional coordinates of atoms',
            nmin=nmin, nmax=nmax)
        poslines = self.caslines[lposn + 3:lposn + 3 + self.Nions]
        for i in range(self.Nions):
            posns[i, :] = [float(p) for p in poslines[i].split()[3:6]]
        return posns

    # Get the external cell pressure used in the calculation
    def get_ext_press(self):
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
        except IndexError:
            press = [0.0]*6
        return press
    
    # Get the Cell Constraints for a Geometry Optimisation
    def get_cell_constrs(self):
        try:
            lconstrs = stri.strindex(self.caslines, 'Cell constraints are:')
            cellconstrs = [int(c) for c in
                           self.caslines[lconstrs].split()[3:9]]
        except UnboundLocalError:
            cellconstrs = None
        return cellconstrs

    def get_cell(self, iteration=-1):
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

    def get_enthalpy(self,iteration=-1):
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract enthalpy.')
        if self.task == 'single':
            raise NotImplementedError(
                'get_enthalpy not implemented if task == single.')
        elif self.task == 'geometry':
            (nmin,nmax)=self.geomrange(iteration=iteration)
            lenthalpy=stri.strindex(self.caslines,'with enthalpy=',
                                    nmin=nmin,nmax=nmax)
            enthalpy=float(self.caslines[lenthalpy].split()[6])
        return enthalpy

    def get_energy(self,iteration=-1):
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract enthalpy.')
        if self.task == 'single':
            (nmin,nmax)=(0,self.Nlines)
        elif self.task == 'geometry':
            (nmin,nmax)=self.geomrange(iteration=iteration)
        try:
            lenergy=stri.strindex(self.caslines,'Final energy, E',
                                  nmin=nmin,nmax=nmax)
            energy=float(self.caslines[lenergy].split()[4])
        except UnboundLocalError:
            lenergy=stri.strindex(self.caslines,'Final energy =',
                                  nmin=nmin,nmax=nmax)
            energy=float(self.caslines[lenergy].split()[3])
        return energy
    
    def get_forces(self,iteration=-1):
        forces=np.zeros((self.Nions,3))
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract forces.')
        if self.task == 'single':
            (nmin,nmax)=(0,self.Nlines)
        elif self.task == 'geometry':
            (nmin,nmax)=self.geomrange(iteration=iteration)
        lforce=stri.eitheror_strsindex(self.caslines,
                                       ['* Forces *','* Symmetrised Forces *']
                                       ,nmin=nmin,nmax=nmax)
        forcelines=self.caslines[lforce+6:lforce+6+self.Nions]
        for i in range(self.Nions):
            newforcelinesplt=[f for f in forcelines[i].split() if f!='(mixed)']
            forces[i,:]=[float(p.replace('(cons\'d)', ''))
                         for p in newforcelinesplt[3:6]]        
        return forces

    def get_stresses(self,iteration=-1):
        stresses=np.zeros((3,3))
        if self.Niterations == 0:
            raise IndexError(
                'No SCF calculations completed so cannot extract stresses.')
        if self.task == 'single':
            (nmin,nmax)=(0,self.Nlines)
        elif self.task == 'geometry':
            (nmin,nmax)=self.geomrange(iteration=iteration)
        lstress=stri.eitheror_strsindex(self.caslines,
                                        ['* Stress Tensor *',
                                         '* Symmetrised Stress Tensor *'],
                                        nmin=nmin,nmax=nmax)
        stresslines=self.caslines[lstress+6:lstress+9]
        for i in range(3):
            stresses[i,:]=[float(p) for p in stresslines[i].split()[2:5]]
        return stresses

    def get_Smax(self,iteration=-1):
        # I know CASTEP computes this information but extracting it was
        # slightly problematic and I'm lazy - it's a very small calculation!
        stresses=self.get_stresses(iteration=iteration)
        abstresses=[]
        for i in range(3):
            S=np.linalg.norm(stresses[i,:])
            abstresses.append(S)
        Smax=max(abstresses)
        return Smax
    
    def get_Fmax(self,iteration=-1):
        """
        if self.task == 'geometry':
            (nmin,nmax)=self.geomrange(iteration=iteration)
            lFmax=stri.strindex(self.caslines,'|F|max',nmin=nmin,nmax=nmax)
            Fmax=float(self.caslines[lFmax].split()[3])
        else:
        """
        # I know CASTEP computes this information but extracting it was
        # slightly problematic and I'm lazy - it's a very small calculation!
        if (1 == 1):
            forces=self.get_forces(iteration=iteration)
            abforces=[]
            for i in range(self.Nions):
                F=np.linalg.norm(forces[i,:])
                abforces.append(F)
            Fmax=max(abforces)
        return Fmax
