import ase
import numpy as np

"""
Module managing conversion between solid solution structures (using the VCA)
and structures with only a single atom per site (e.g. for phonon calculations).

Two formats for mixkeys may be considered:

"site" mixkey (sites in fractional coordinates) e.g.
sitemixkey = {'0.5 0.5 0.5': ('Ca', {'Ca': 0.5, 'Sr': 0.5}),
              '0.0 0.0 0.0': ('Ge', {'Ge': 1.0}),
              '0.5 0.5 0.0': ('O', {'O': 1.0}),
              '0.0 0.5 0.5': ('O', {'O': 1.0}),
              '0.5 0.0 0.5': ('O', {'O': 1.0})}

"element" mixkey, e.g.
elemmixkey = {'Ca': {'Ca': 0.5, 'Sr': 0.5},
              'Ge': {'Ge': 1.0},
              'O': {'O': 1.0}}
WARNING: be careful with labelling if this format is used when different
sites exist containing the same element -- key labels in the pure structure
must be unique
"""

# Defaults -- these are relevant when printing cell files
# Note that these are sensible defaults for A2BO4 Ruddlesden-Popper oxides
pseudo_default = True
kpts_default = [8, 8, 4]
kpts_offset_default = [0.1, 0.1, 0.5]
spins_default = None
pressure_default = [0.0]*6
cell_constrs_default = [1, 2, 3, 0, 0, 0]  # Assumes orthorhombic cell


def create_mixture(atoms, elem_mixkey):
    """
    ase.Atoms atoms : atomic structure with no mixing (one atom per site)
    mixkey elem_mixkey : mapping of the elem_mixkey format
    
    returns
    ase.Atoms mixatoms : atomic structure with multiple atoms per site """
    posns = atoms.get_positions()
    elems = atoms.get_chemical_symbols()
    newposns = posns.copy()
    newelems = [] + elems
    for key in list(elem_mixkey.keys()):
        subs = [elem for elem in elem_mixkey[key] if elem != key]
        if len(subs):
            keyinds = [i for i, elem in enumerate(elems) if elem == key]
            keyposns = posns[keyinds]
            for sub in subs:
                newposns = np.vstack((newposns, keyposns))
                newelems += [sub]*len(keyinds)
    cell = atoms.get_cell()
    mixatoms = ase.Atoms(symbols=newelems, positions=newposns,
                         cell=cell, pbc=True)
    return mixatoms


class mixmap():
    """ Class used for mapping between structures with mixed atoms and pure
    atoms on a single site."""
    
    def __init__(self, mixatoms, mixkey, wttol=0.0001, postol=0.5):
        """ should be initialised for a particular mixed atom structure
        
        ase.Atoms mixatoms : atomic structure with multiple atoms on same site
        dict mixkey : info atom mix per site (can be site or elem format)
        float wttol : weights should sum to 1.0 (tolerance for rounding errors)
        float postol : look for atomic site keys within this tolerance of posn
        """
        self.check_wts(mixkey, wttol)  # check that site weights sum to 1.0
        self.mixkey = mixkey
        self.postol = postol
        # This info is fixed for this mixmap instance
        self.mixelems = mixatoms.get_chemical_symbols()
        self.mixions = mixatoms.get_number_of_atoms()
        self.mixmasses = mixatoms.get_masses()
        
        # This sets up mappings between pure and mix structures
        # And pure structure info (elements, Nions, masses)
        self.setup_maps(mixatoms)
        
        self.setcellparams()  # Initialise calc. params at default values
    
    def setup_maps(self, mixatoms):
        """ Setup pure2mix_map and mix2pure_map mappings.
        
        These mappings map an atom index in the pure structure to indices in
        the mix structure or visa versa.
        This method also sets up attributes of the pure structure that
        remain fixed for this instance (i.e. Nions, elems, masses)
        
        ase.Atoms mixatoms : structure with multiple atoms on same site"""
        
        mixposns = mixatoms.get_positions()
        cell = mixatoms.get_cell()
        
        pureelems = []     # Element list for pure structure
        puremasses = []    # Masses for each pure site (average of mix atoms)
        
        pure2mix_map = {}  # Map pure indices to mix indices
        mix2pure_map = {}  # Map mix indices to a pure index
        mixsitemixes = {}  # Gives the weight corresponding to each mix index
        
        p = 0  # Index of atoms in the pure structure
        m = 1  # Index of mixture atoms
        for i in range(self.mixions):
            mixelem = self.mixelems[i]
            mixposn = mixposns[i, :]
            matchkey = self.sitematch(mixelem, mixposn, cell)
            if len(matchkey.split()) == 3:
                pureelem, wts = self.mixkey[matchkey]
            else:
                pureelem = matchkey
                wts = self.mixkey[matchkey]
            if len(wts) == 1:  # This is the trivial case
                puremasses += [self.mixmasses[i]]
                pureelems += [pureelem]
                mix2pure_map[i] = p
                pure2mix_map[p] = {mixelem: (1.0, i)}
                mixsitemixes[i] = (0, 1.0)
                p += 1
            
            elif mixelem == pureelem:  # Mixed atoms on this site
                mass = 0
                sitemapdict = {}
                # Cycle through all mixture elements on site
                for elem in list(wts.keys()):
                    for j in range(self.mixions):
                        posdiff = np.linalg.norm(mixposns[j]-mixposn)
                        if self.mixelems[j] == elem and posdiff < self.postol:
                            mass += self.mixmasses[j]*wts[elem]
                            sitemapdict[elem] = (wts[elem], j)
                            mixsitemixes[j] = (m, wts[elem])
                puremasses.append(mass)
                pureelems.append(mixelem)
                mix2pure_map[i] = p
                pure2mix_map[p] = sitemapdict
                p += 1
                m += 1
        
        self.pureions = len(pureelems)
        self.pureelems = pureelems
        self.puremasses = np.array(puremasses)
        self.pure2mix_map = pure2mix_map
        self.mix2pure_map = mix2pure_map
        self.mixsitemixes = mixsitemixes
    
    @staticmethod
    def check_site_mixkey(site_mixkey):
        """ Raises an error if the input mixkey is not of the site format """
        error = True
        if isinstance(site_mixkey, dict):
            layer1 = list(site_mixkey.values())[0]
            cont = False
            if isinstance(layer1, (tuple, list)):
                layer2 = layer1[1]
                cont = True
            elif isinstance(layer1, dict):
                # This could still be unambiguous
                layer2 = list(layer1.values())[0]
                cont = True
            if cont:
                if isinstance(layer2, dict):
                    layer3 = list(layer2.values())[0]
                    if isinstance(layer3, float):
                        error = False
        if error:
            raise ValueError("mixkey not of site form: " +
                             "{'x y z': ('A': {'A': 0.4, 'B':0.6, ...}), ...}")
    
    @staticmethod
    def check_elem_mixkey(elem_mixkey):
        """ Raises an error if the input mixkey is not of the elem format """
        error = True
        if isinstance(elem_mixkey, dict):
            layer1 = list(elem_mixkey.values())[0]
            if isinstance(layer1, dict):
                layer2 = list(layer1.values())[0]
                if isinstance(layer2, float):
                    error = False
        if error:
            raise ValueError("mixkey not of elem form: " +
                             "{'A': {'A': 0.4, 'B':0.6, ...}, ...}")
    
    @staticmethod
    def site2elem_mixkey(site_mixkey):
        """ Converts a mixkey of the site format to one of the elem format """
        mixmap.check_site_mixkey(site_mixkey)
        elem_mixkey = {}
        for sitekey in site_mixkey:
            siteelem, sitedict = site_mixkey[sitekey]
            if siteelem not in elem_mixkey:
                elem_mixkey[siteelem] = sitedict
            else:
                if sitedict != elem_mixkey[siteelem]:
                    raise ValueError("Cannot convert site_mixkey to " +
                                     "elem_mixkey because different mixtures "
                                     + "for same  element on different sites.")
        return elem_mixkey
    
    @staticmethod
    def closestimage(a, b, lats, rtn_dist=False):
        """ Give the coordinates of the closest atom b to atom a given periodic
        boundary conditions
        np.array(3) a : absolute coordinates of reference atom (to be close to)
        np.array(3) b : absolute coordinates of movable atom
        np.array(3, 3) lats : unit cell vector
        bool dist : return a tuple of (bprime, dist) not just bprime
        
        returns
        np.array(3) bprime : absolute coordinates of closest image of b to a
        """
        dists = []
        bprimes = np.zeros((27, 3))
        x = 0
        for h in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                for l in [-1, 0, 1]:
                    bprime = b + h*lats[0, :] + k*lats[1, :] + l*lats[2, :]
                    bprimes[x, :] = bprime
                    dists.append(np.linalg.norm(a-bprime))
                    x += 1
        dist = min(dists)
        bprime = bprimes[dists.index(dist), :]
        if rtn_dist:
            return bprime, dist
        else:
            return bprime
    
    @staticmethod
    def check_wts(mixkey, wttol=0.0001):
        """ Check that the atom weights for each site sums to 1.0 """
        for sitekey in list(mixkey.keys()):
            if len(sitekey.split()) == 3:  # site_mixkey
                wts = mixkey[sitekey][1]
            elif len(sitekey.split()) == 1:  # elem_mixkey
                wts = mixkey[sitekey]
            else:
                raise KeyError(sitekey+' not recognised as mixkey.')
            sitewt = sum(list(wts.values()))
            if abs(sitewt-1) > wttol:
                raise AttributeError('Sum of concs on site ' + sitekey
                                     + ' equals ' + str(sitewt) +
                                     ' which does not make sense.')
    
    def sitematch(self, elem, posn, cell):
        """
        Match a mix atomic site to that in the mixkey either by site or
        element (depending on format of mixkey)
        
        str elem : element name
        np.array(3) posn : absolute position of atom
        np.array(3, 3) cell : unit cell vectors
        
        returns
        str matchkey : key from self.mixkeys that matches element and site
        """
        sitekeys = list(self.mixkey.keys())
        matchkey = None
        if len(sitekeys[0].split()) == 3:  # site format for mixkey
            dists = []
            for sitekey in sitekeys:
                site = np.array([float(s) for s in sitekey.split()])
                siteposn = np.dot(site, cell)
                vector, dist = self.closestimage(posn, siteposn, cell,
                                                 rtn_dist=True)
                dists += [dist]
            matchkey = sitekeys[dists.index(min(dists))]
        elif len(sitekeys[0].split()) == 1:  # elem format for mixkey
            for sitekey in sitekeys:
                if elem in self.mixkey[sitekey]:
                    matchkey = sitekey
            if matchkey is None:
                raise KeyError('Element: '+elem+' does not appear in mixkeys')
        else:
            raise KeyError(str(sitekeys[0]) + '  not recognised as mixkey.')
        return matchkey
    
    def pure2mix(self, pureatoms):
        """ Convert a pure ase.Atoms structure to a mixed structure """
        pureposns = pureatoms.get_positions()
        cell = pureatoms.get_cell()
        mixposns = np.zeros((self.mixions, 3))
        for i in range(self.pureions):
            sites = self.pure2mix_map[i]
            for siteelem in list(sites.keys()):
                mixposns[sites[siteelem][1], :] = pureposns[i, :]
        mixatoms = ase.Atoms(positions=mixposns, symbols=self.mixelems,
                             cell=cell, pbc=True)
        return mixatoms
    
    def mix2pure(self, mixatoms, phonopy=False):
        """ Convert a mix ase.atoms structure to a pure structure
        bool phonopy : if True will return a phonopy atoms object not ase """
        mixposns = mixatoms.get_positions()
        cell = mixatoms.get_cell()
        pureposns = np.zeros((self.pureions, 3))
        for i in range(self.mixions):
            try:
                pureposns[self.mix2pure_map[i], :] = mixposns[i, :]
            except KeyError:
                pass
        if phonopy:
            from phonopy.structure import atoms
            pureatoms = atoms.PhonopyAtoms(positions=pureposns,
                                           symbols=self.pureelems, cell=cell,
                                           masses=self.puremasses, pbc=True)
        else:
            pureatoms = ase.Atoms(positions=pureposns, symbols=self.pureelems,
                                  cell=cell, masses=self.puremasses, pbc=True)
        return pureatoms
    
    def pinch_posns(self, mixatoms):
        """ Ensure that all mix atoms that are meant to occupy the same
        site actually have the same coordinate -- can be an issue in
        CASTEP geometry relaxation if structure is polar. """
        mixposns = mixatoms.get_positions()
        cell = mixatoms.get_cell()
        pureposns = np.zeros((self.pureions, 3))
        for i in range(self.pureions):
            sitedict = self.pure2mix_map[i]
            posn = np.zeros((3))
            for key in list(sitedict.keys()):
                wt, idx = sitedict[key]
                posn += self.closestimage(posn, mixposns[idx], cell)*wt
            pureposns[i, :] = posn
        pureatoms = ase.Atoms(positions=pureposns,
                              symbols=self.pureelems,
                              cell=cell, masses=self.puremasses,
                              pbc=True)
        pureatoms.wrap()
        return self.pure2mix(pureatoms)
    
    def mix2pure_spins(self, mixspins):
        """ Convert mixed spins to pure spins (not accounting for weights!) """
        purespins = np.zeros((self.pureions))
        for i in range(self.mixions):
            try:
                purespins[self.mix2pure_map[i]] = mixspins[i]
            except KeyError:
                pass
        return purespins
    
    def pure2mix_spins(self, purespins):
        """ Convert pure spins to mixed spins (assumes they are equal) """
        mixspins = np.zeros((self.mixions))
        for i in range(self.pureions):
            sites = self.pure2mix_map[i]
            for siteelem in list(sites.keys()):
                mixspins[sites[siteelem][1]] = purespins[i]
        return mixspins
    
    def mix2pure_forces(self, mixforces):
        """ Convert forces from mixed structure to pure forces (by taking
        a weighted average on each site) """
        pureforces = np.zeros((self.pureions, 3))
        for i in range(self.mixions):
            try:
                pureforces[self.mix2pure_map[i], :] = mixforces[i, :]
            except KeyError:
                pass
        return pureforces
    
    def casprint(self, atoms, cellfile, pure=False):
        """ Produce a CASTEP .cell file from an atoms object
        
        ase.Atoms atoms : structure to produce file from
        str cellfile : filename (incl. path) to write to
        bool pure : if True, produce .cell of pure struc, otherwise mix """
        
        # Print cell
        caslines = ['%BLOCK LATTICE_CART\n']
        cell = atoms.get_cell()
        for i in range(3):
            caslines += ['\t'+'\t'.join([str('{0:.8f}'.format(cell[i, j]))
                                         for j in range(3)])+'\n']
        caslines += ['%ENDBLOCK LATTICE_CART\n']+['\n']
        
        # Print cell constraints (if not [1, 2, ... 6])
        if any([self.cell_constrs[i] != i+1 for i in range(6)]):
            caslines += ['%BLOCK cell_constraints\n']
            caslines += ['\t'+'\t'.join([str(int(self.cell_constrs[j]))
                                         for j in range(3)])+'\n']
            caslines += ['\t'+'\t'.join([str(int(self.cell_constrs[j]))
                                         for j in range(3, 6)])+'\n']
            caslines += ['%ENDBLOCK cell_constraints\n']+['\n']
                
        # Print elements, atomic positions, spins and mix weights
        mixelems = atoms.get_chemical_symbols()
        Natoms = len(mixelems)
        spins = self.spins
        if pure:
            spins = self.mix2pure_spins(spins)
            Natoms = self.pureions
        if self.frac:
            caslines += ['%BLOCK POSITIONS_FRAC\n']
            mixposns = atoms.get_scaled_positions()
        else:
            caslines += ['%BLOCK POSITIONS_ABS\n']
            mixposns = atoms.get_positions()
        for i in range(Natoms):
            (m, wt) = self.mixsitemixes[i]
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
            spin = spins[i]
            if spin == 0:
                spinstring = ''
            else:
                spinstring = '\tSPIN='+str(spin)
            posstring = '\t'.join(
                [str('{0:.8f}'.format(mixposns[i, j])) for j in range(3)])
            caslines += ['\t' + mixelems[i] + '\t' + posstring
                         + spinstring + mixstring + '\n']
        if self.frac:
            caslines += ['%ENDBLOCK POSITIONS_FRAC\n']
        else:
            caslines += ['%ENDBLOCK POSITIONS_ABS\n']

        # Print k-points and offset
        caslines += ['\n']
        caslines += ['kpoints_mp_grid = '+' '.join([str(self.kpoints[j])
                                                    for j in range(3)])+'\n']
        caslines += ['kpoints_mp_offset = '+' '.join(
            [str(self.kpoints_offset[j]) for j in range(3)])+'\n']
        
        # Pseudo potentials
        caslines += ['\n']+['%BLOCK SPECIES_POT\n']
        for elem in list(self.pseudos.keys()):
            caslines += ['\t'+elem+' '+self.pseudos[elem]+'\n']
        caslines += ['%ENDBLOCK SPECIES_POT\n']+['\n']
        
        # Symmetry statements
        if self.sym_gen:
            caslines += ['symmetry_generate\n']+['\n']
        if self.snap_sym:
            caslines += ['snap_to_symmetry\n']+['\n']
        
        # External pressure
        if any([self.pressure[i] != 0.0 for i in range(6)]):
            caslines += ['%BLOCK external_pressure\n']+['\tGPA\n']
            caslines += ['\t'+'\t'.join(
                [str('{0:.8f}'.format(self.pressure[j]))
                 for j in [0, 5, 4]])+'\n']
            caslines += ['\t\t\t'+'\t'.join([str('{0:.8f}'.format(
                self.pressure[j])) for j in [1, 3]])+'\n']
            caslines += ['\t\t\t\t\t'+'\t'.join([str('{0:.8f}'.format(
                self.pressure[j])) for j in [2]])+'\n']
            caslines += ['%ENDBLOCK external_pressure\n']+['\n']
        
        # Ionic constraints
        # WARNING: at the moment this doesn't switch pureatoms <--> mixatoms
        if self.ion_constrs is not None:
            caslines += ['%BLOCK IONIC_CONSTRAINTS\n']
            k = 1
            for i, elem in enumerate(mixelems):
                for j in range(3):
                    zeros = [0.0, 0.0, 0.0]
                    if self.ion_constrs[i, j]:
                        zeros[j] = 1.0
                        caslines += [str(k) + '\t' + elem + '\t'
                                     + str(i+1) + '\t' +
                                     '\t'.join([str(z) for z in zeros])+'\n']
                        k += 1
            caslines += ['%ENDBLOCK IONIC_CONSTRAINTS\n']+['\n']
        
        # This writes the .cell file
        open(cellfile, 'w').writelines(caslines)
    
    def setcellparams(self, pseudos=pseudo_default, frac=True,
                      kpoints=kpts_default, kpoints_offset=kpts_offset_default,
                      sym_gen=True, snap_sym=True, spins=spins_default,
                      pressure=pressure_default,
                      cell_constrs=cell_constrs_default, ion_constrs=None):
        """ Set additonal info for the CASTEP .cell file """
        if pseudos:
            elemset = list(set(self.mixelems))
            pseudos = {"Mg": "1|1.8|3|4|4|30N:31L:32N",
                    "O": "2|1.2|23|26|31|20NN:21NN(qc=9)",
                    "La": "4|2.1|15|17|20|50N:60N:51N:52N:43N{4f0.1}(qc=6,q3=8)",
                    "Al": "1|1.6|6|7|8|30N:31L:32N",
                    "Ba":"2|2.0|8|10|11|50N:60N:51N(qc=5.5)"}
        self.pseudos = pseudos
        self.frac = frac
        self.kpoints = kpoints
        self.kpoints_offset = kpoints_offset
        self.sym_gen = sym_gen
        self.snap_sym = snap_sym
        if spins is None:
            spins = [0]*self.mixions
        self.spins = spins
        self.pressure = pressure
        self.cell_constrs = cell_constrs
        self.ion_constrs = ion_constrs
