import ase
import numpy as np
from phonopy.structure import atoms # PHONOPY unfortunately cannot take regular Atoms objects
import string

# Class for performing phonon calculations (using PHONOPY) on CASTEP strucutres with mixed atoms (within the VCA)

# Defaults
pseudo_default=True
kpts_default=[5,5,1]
kpts_offset_default=[0.1,0.1,0.5]
spins_default=None
wtdefault=0.0001   # Each site should have an exact weight of 1.0, this tol is just to avoid numerical rounding errors 
posdefault=0.5  # It is STRONGLY advised to setup the mixmap instance with a fully consistent geometry, once you have the mapping you can resuse later with no troubles
sitemixkeydefault={'A':{'Ca':{'Ca':1.0}},'B':{'Ge':{'Ge':1.0}},'O':{'O':{'O':1.0}}}  # Never rely on this since clearly there is no mixing (no longer there is default)
pressure_default=[0.0]*6
cell_constrs_default=[1,1,3,0,0,0]  # Note this default is for a tetragonal cell!

def mixkey2sitemixkey(mixkey):
    i = 0
    sitemixkey = {}
    knownelems = []
    for sitekey in mixkey:
        elems = sorted(mixkey[sitekey].keys())
        if any([elem in knownelems for elem in elems]):
            try:
                for letter in list(sitemixkey.keys()):
                    if any([key in elems for key in sitemixkey[letter].keys()]):
                        for key in list(sitemixkey[letter].keys()):
                            if key in elems:
                                if ((len(elems) == 1 and not
                                     sitemixkey[letter][key] == mixkey[sitekey][elems[0]]) or
                                    (len(elems) > 1 and not sitemixkey[letter][key] == mixkey[sitekey])): 
                                    raise ValueError('Mix weights do not match: '
                                                     +str(sitemixkey[letter][key])+'\t'+str(mixkey[sitekey]))
            except:
                raise ValueError('Mixkey cannot be converted to sitemixkey.')
        else:
            knownelems += elems
            if len(elems) == 1:
                sitemixkey[string.ascii_uppercase[i]] = mixkey[sitekey]
            else:
                sitemixkey[string.ascii_uppercase[i]] = {elems[0]:mixkey[sitekey]}
            i += 1
    return sitemixkey

def mix_atoms(atoms, sitemixkey):
    Nats = len(atoms)
    posns = atoms.get_positions()
    elems = atoms.get_chemical_symbols()
    newposns = posns.copy()
    newelems = [] + elems
    for letter in list(sitemixkey.keys()):
        for key in list(sitemixkey[letter].keys()):
            subs = [elem for elem in sitemixkey[letter][key] if elem != key]
            if len(subs):
                keyinds = [i for i, elem in enumerate(elems) if elem == key]
                keyposns = posns[keyinds]
                for sub in subs:
                    newposns = np.vstack((newposns, keyposns))
                    newelems += [sub]*len(keyinds)
    cell = atoms.get_cell()
    casatoms = ase.Atoms(symbols=newelems, positions=newposns,
                         cell=cell, pbc=True)
    return casatoms

def closestimage(a,b,lats):
    dists=[]
    bprimes=np.zeros((27,3))
    x=0
    for h in [-1,0,1]:
        for k in [-1,0,1]:
            for l in [-1,0,1]:
                bprime=b+h*lats[0,:]+k*lats[1,:]+l*lats[2,:]
                bprimes[x,:]=bprime
                dists.append(np.linalg.norm(a-bprime))
                x+=1
    return bprimes[dists.index(min(dists)),:]

def closestdist(a,b,lats):
    dists=[]
    for h in [-1,0,1]:
        for k in [-1,0,1]:
            for l in [-1,0,1]:
                bprime=b+h*lats[0,:]+k*lats[1,:]+l*lats[2,:]
                dists.append(np.linalg.norm(a-bprime))
    return min(dists)

# Class instance is really a recipe for mapping between ase.Atoms objects with (cas) and without (phon) mix atoms
class mixmap():
    # Should be initialised for a particular cas structure
    def __init__(self, casatoms, mixkey, wttol=wtdefault, postol=posdefault):
        # mixkey
        self.mixkey=mixkey
        self.postol=postol
        # Check that the mix for each site sums to 1.0
        for sitekey in list(mixkey.keys()):
            #print(mixkey[sitekey].values()) 
            wts=list(mixkey[sitekey].values())[0]
            #print wts
            sitewt=sum(list(wts.values()))
            if sitewt<1.0-wttol or sitewt>1.0+wttol:
                raise AttributeError('Sum of concs on site '+sitekey+' equals '+str(sitewt)+' which does not make sense.')
        self.caselems=casatoms.get_chemical_symbols()
        self.casions=casatoms.get_number_of_atoms()
        casposns=casatoms.get_positions() # no need to make this an attribute since not used outside initialisation
        casmasses=casatoms.get_masses() # no need to make this an attribute since not used outside initialisation
        cell=casatoms.get_cell()
        cassitemixes={}
        p2cmap={}
        c2pmap={}
        phonelems=[]
        phonposns=[]
        phonmasses=[]
        p=0  # Index to keep track of the phon mapping
        m=1  # Index to keep track of the mixture atoms (start counting at 1)
        for i in range(self.casions):  # i is the index of the cas ion remember
            caselem=self.caselems[i]
            casposn=casposns[i,:]
            matchkey=self.sitematch(caselem,casposn,cell)
            phonelem=list(self.mixkey[matchkey].keys())[0]
            #print(str(i)+'\t'+caselem+'\t'+phonelem+'\t'+str(casposn)+'\t'+
            #      str(matchkey))
            wts=self.mixkey[matchkey][phonelem]
            if len(wts.keys())>1: # if pure element no need for mix (actually doesn't work with mix!!)
                if caselem==phonelem:   # ignore all mix elements in the outer loop not to be included in phon 
                    mass=0
                    sitemapdict={}
                    # Cycle through all secondary mixture elements for this site
                    for secondelem in list(wts.keys()):
                        for j in range(self.casions):
                            if self.caselems[j]==secondelem and np.linalg.norm(casposns[j]-casposn)<postol:
                                mass+=casmasses[j]*wts[secondelem]
                                sitemapdict[secondelem]=(wts[secondelem],j)
                                cassitemixes[j]=(m,wts[secondelem])
                    phonmasses.append(mass)
                    phonelems.append(caselem)
                    phonposns.append(casposn)
                    c2pmap[i]=p
                    p2cmap[p]=sitemapdict
                    p+=1
                    m+=1
            else:
                phonmasses.append(casmasses[i])
                phonelems.append(phonelem)
                phonposns.append(casposn)
                c2pmap[i]=p
                p2cmap[p]={caselem:(1.0,i)}
                cassitemixes[i]=(0,1.0)
                p+=1
        self.phonions=len(phonelems)
        self.phonelems=phonelems
        self.phonmasses=np.array(phonmasses)
        self.phonposns=np.array(phonposns)
        self.p2cmap=p2cmap
        self.c2pmap=c2pmap
        self.cassitemixes=cassitemixes
        self.setcascellinfo() # Set the default info for printing cell files (can be overwritten later)

    
    # Match a cas atomic site to that in the mixkey either by posn or element (depending on format of key)
    def sitematch(self, elem, posn, cell):
        sitekeys=list(self.mixkey.keys())
        matchkey=''
        # Can be atom by atom: atommixkey
        # e.g. atommixkey = {'0.0 0.0 0.0':{'Ge':{'Ge':0.8,'Sr':0.2}}, ... }
        if len(sitekeys[0].split())==3:
            dists = []
            for sitekey in sitekeys:
                siteposn=np.dot(np.array([float(s) for s in sitekey.split()]),cell)
                #dists += [np.linalg.norm(posn-siteposn)]
                dists += [closestdist(posn,siteposn,cell)]
                # This isn't technically correct since siteposn labels are in frac
                #if np.linalg.norm(posn-siteposn)<self.postol:
                #    matchkey=sitekey
            matchkey = sitekeys[dists.index(min(dists))]
            #if matchkey=='':
            #    raise KeyError('Position: ['+' '.join([str(posn[j]) for j in range(3)])+'] does not appear in mixkeys')
        # Can be site specific: sitemixkey
        # e.g. sitemixkey = {'B':{'Ge':{'Ge':0.8,'Sr':0.2}}, ... }
        elif len(sitekeys[0].split())==1:
            for sitekey in sitekeys:
                if elem in list(list(self.mixkey[sitekey].values())[0].keys()):   # elem appears as any of the elements that mix on that site
                    matchkey=sitekey
            if matchkey=='':
                raise KeyError('Element: '+elem+' does not appear in mixkeys')
        else:
            raise TypeError('Type: '+type(sitekeys[0])+' of mixkey variable not recognised as mixkey.')
        return matchkey
            
    # Convert a PHONOPY atoms object to a CASTEP atoms object
    def phon2cas(self, phonatoms):
        phonposns=phonatoms.get_positions()
        cell=phonatoms.get_cell()
        casposns=np.zeros((self.casions,3))
        for i in range(self.phonions):
            sites=self.p2cmap[i]
            for siteelem in list(sites.keys()):
                casposns[sites[siteelem][1],:]=phonposns[i,:]
        casatoms=ase.Atoms(positions=casposns,symbols=self.caselems,cell=cell,pbc=True)        
        return casatoms
    
    # Convert a CASTEP atoms object to an ase/PHONOPY atoms object
    def cas2phon(self, casatoms, phonopy=False):
        casposns=casatoms.get_positions()
        cell=casatoms.get_cell()
        phonposns=np.zeros((self.phonions,3))
        for i in range(self.casions):
            try:
                phonposns[self.c2pmap[i],:]=casposns[i,:]
            except KeyError:
                pass
        if phonopy:
            phonatoms = atoms.PhonopyAtoms(positions=phonposns,
                                           symbols=self.phonelems,
                                           cell=cell, masses=self.phonmasses,
                                           pbc=True)
        else:
            phonatoms = ase.Atoms(positions=phonposns,
                                  symbols=self.phonelems,
                                  cell=cell, masses=self.phonmasses,
                                  pbc=True)            
        return phonatoms

    def pinch_posns(self, casatoms):
        casposns = casatoms.get_positions()
        cell = casatoms.get_cell()
        phonposns = np.zeros((self.phonions,3))
        for i in range(self.phonions):
            sitedict = self.p2cmap[i]
            posn = np.zeros((3))
            for key in list(sitedict.keys()):
                wt, idx = sitedict[key]
                posn += closestimage(posn, casposns[idx], cell)*wt
            phonposns[i, :] = posn        
        phonatoms = ase.Atoms(positions=phonposns,
                              symbols=self.phonelems,
                              cell=cell, masses=self.phonmasses,
                              pbc=True)
        phonatoms.wrap()
        return self.phon2cas(phonatoms)
    
    def cas2phon_spins(self, casspins):
        phonspins = np.zeros((self.phonions))
        for i in range(self.casions):
            try:
                phonspins[self.c2pmap[i]]=casspins[i]
            except KeyError:
                pass
        return phonspins

    def phon2cas_spins(self, phonspins):
        casspins = np.zeros((self.casions))
        for i in range(self.phonions):
            sites = self.p2cmap[i]
            for siteelem in list(sites.keys()):
                casspins[sites[siteelem][1]] = phonspins[i]
        """
        for i in range(self.casions):
            casspins[i] = phonspins[self.c2pmap[i]]
        """
        return casspins

    # Convert a CASTEP atoms object to a PHONOPY atoms object
    def cas2phon_forces(self, casforces):
        phonforces=np.zeros((self.phonions,3))
        for i in range(self.casions):
            try:
                phonforces[self.c2pmap[i],:]=casforces[i,:]
            except KeyError:
                pass
        return phonforces

    # Print a CASTEP .cell file from a cas atoms object
    def casprint(self, atoms, cellfile, phon=False):
        caslines=['%BLOCK LATTICE_CART\n']
        cell=atoms.get_cell()
        for i in range(3):
            caslines+=['\t'+'\t'.join([str('{0:.8f}'.format(cell[i,j]))
                                       for j in range(3)])+'\n']
        caslines+=['%ENDBLOCK LATTICE_CART\n']+['\n']
        # if numbered 1 to 6, vary cell => default behaviour
        if any([self.cell_constrs[i]!=i+1 for i in range(6)]):
            caslines+=['%BLOCK cell_constraints\n']
            caslines+=['\t'+'\t'.join([str(int(self.cell_constrs[j]))
                                       for j in range(3)])+'\n']
            caslines+=['\t'+'\t'.join([str(int(self.cell_constrs[j]))
                                       for j in range(3,6)])+'\n']
            caslines+=['%ENDBLOCK cell_constraints\n']+['\n']
        caselems=atoms.get_chemical_symbols()
        Natoms=len(caselems)
        spins=self.spins
        if phon:
            spins=self.cas2phon_spins(spins)
        if self.frac==1:
            caslines+=['%BLOCK POSITIONS_FRAC\n']
            casposns=atoms.get_scaled_positions()
            for i in range(Natoms):
                (m,wt)=self.cassitemixes[i]
                if phon:
                    wt=1.0
                if wt==1.0:
                    mixstring=''
                elif wt<=0.0 or wt>1.0:
                    raise ValueError('Trying to print ion index '+str(i)+
                                     ' with weight '+str(wt)+
                                     ' and do not know what to do.')
                else:
                    mixstring='\tMIXTURE=('+str(m)+' '+str(wt)+')'
                spin=spins[i]
                if spin==0:
                    spinstring=''
                else:
                    spinstring='\tSPIN='+str(spin)
                caslines+=['\t'+caselems[i]+'\t'+'\t'.join([str('{0:.8f}'.format(casposns[i,j])) for j in range(3)])+spinstring+mixstring+'\n']
            caslines+=['%ENDBLOCK POSITIONS_FRAC\n']
        else:
            caslines+=['%BLOCK POSITIONS_ABS\n']
            casposns=atoms.get_positions()
            for i in range(Natoms):
                (m,wt)=self.cassitemixes[i]
                if wt==1.0:
                    mixstring=''
                elif wt<=0.0 or wt>1.0:
                    raise ValueError('Trying to print ion index '+str(i)+' with weight '+str(wt)+' and do not know what to do.')
                else:
                    mixstring='\tMIXTURE=('+str(m)+' '+str(wt)+')'
                spin=spins[i]
                if spin==0:
                    spinstring=''
                else:
                    spinstring='\tSPIN='+str(spin)
                caslines+=['\t'+caselems[i]+'\t'+'\t'.join([str('{0:.8f}'.format(casposns[i,j])) for j in range(3)])+spinstring+mixstring+'\n']
            caslines+=['%ENDBLOCK POSITIONS_ABS\n']
        caslines+=['\n']
        caslines+=['kpoints_mp_grid = '+' '.join([str(self.kpoints[j])
                                                  for j in range(3)])+'\n']
        caslines+=['kpoints_mp_offset = '+' '.join(
            [str(self.kpoints_offset[j]) for j in range(3)])+'\n']
        caslines+=['\n']+['%block species_pot\n']
        for elem in list(self.pseudos.keys()):
            caslines+=['\t'+elem+' '+self.pseudos[elem]+'\n']
        caslines+=['%endblock species_pot\n']+['\n']
        if self.sym_gen==True:
            caslines+=['symmetry_generate\n']+['\n']
        if self.snap_sym==True:
            caslines+=['snap_to_symmetry\n']+['\n']
        if any([self.pressure[i]!=0.0 for i in range(6)]):
            caslines+=['%BLOCK external_pressure\n']+['\tGPA\n']
            caslines+=['\t'+'\t'.join([str('{0:.8f}'.format(self.pressure[j]))
                                       for j in [0,5,4]])+'\n']
            caslines+=['\t\t\t'+'\t'.join([str('{0:.8f}'.format(
                self.pressure[j])) for j in [1,3]])+'\n']
            caslines+=['\t\t\t\t\t'+'\t'.join([str('{0:.8f}'.format(
                self.pressure[j])) for j in [2]])+'\n']
            caslines+=['%ENDBLOCK external_pressure\n']+['\n']
        # WARNING: at the moment this doesn't switch phonatoms <--> casatoms
        if self.ion_constrs is not None:
            caslines+=['%BLOCK IONIC_CONSTRAINTS\n']
            k = 1
            for i, elem in enumerate(caselems):
                for j in range(3):
                    zeros = [0.0, 0.0, 0.0]
                    if self.ion_constrs[i, j]:
                        zeros[j] = 1.0
                        caslines+=[str(k)+'\t'+elem+'\t'+str(i+1)+'\t'+
                                   '\t'.join([str(z) for z in zeros])+'\n']
                        k += 1
            caslines+=['%ENDBLOCK IONIC_CONSTRAINTS\n']+['\n']
        open(cellfile,'w').writelines(caslines) # Should now have written the castep cell file
            
    # Set additonal info for the CASTEP .cell file
    def setcascellinfo(self, pseudos=pseudo_default, frac=True,
                       kpoints=kpts_default, kpoints_offset=kpts_offset_default,
                       sym_gen=True, snap_sym=True, spins=spins_default,
                       pressure=pressure_default,
                       cell_constrs=cell_constrs_default, ion_constrs=None):
        if pseudos==True: # i.e. default behaviour
            elemset=list(set(self.caselems))
            pseudos={}
            for elem in elemset:
                pseudos[elem]=elem+'_OTF_E1650_PBEsol_NC.usp'
        self.pseudos=pseudos
        self.frac=frac
        self.kpoints=kpoints
        self.kpoints_offset=kpoints_offset
        self.sym_gen=sym_gen
        self.snap_sym=snap_sym
        if spins is None:
            spins=[0]*self.casions
        self.spins=spins
        self.pressure=pressure
        self.cell_constrs=cell_constrs
        self.ion_constrs=ion_constrs
