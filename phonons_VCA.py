import phonopy
import mixmap
import readmixcastep as rc
import numpy as np
import glob
import pinchposns as pinpos

""" Wrapper to manage phonon calculations using phonopy and CASTEP for
solid solutions made with the virtual crystal approximation (VCA) """


def gen_perturbations(casfile, supercell='Gamma'):
    """ Genterate input .cell files for CASTEP singlepoint energy simulations
    as part of a phonopy phonon calculation of a solid solution.
    
    str casfile : filename (incl. path) to .castep file to compute phonons of
    np.array(3, 3) supercell : size of real space supercell dictates k-points
    If supercell is 'Gamma', only unit cell displacements are created """
    # Load the mixed data
    chem = casfile.replace('.castep', '')
    cas = rc.readcas(casfile)
    mixatoms = cas.extract_struc()
    spins = cas.get_final_spin()
    press = cas.get_ext_press()
    
    # Create pure cell
    mixkey = cas.get_mixkey()
    mapping = mixmap.mixmap(mixatoms, mixkey)
    pureatoms = mapping.mix2pure(mixatoms)
    mapping.setcellparams(spins=spins, pressure=press)

    # Generate perturbed cells (displacements) of pure cell
    if supercell == 'Gamma':  # Gamma-point means only single unit cell
        supercell = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    GammaPhonon = phonopy.Phonopy(pureatoms, supercell, symprec=1e-4,
                                  factor=phonopy.units.VaspToCm)
    GammaPhonon.generate_displacements(distance=0.05)
    displcells = GammaPhonon.get_supercells_with_displacements()
    
    # Convert displaced pure cells to mixtures and write files
    for d, displcell in enumerate(displcells):
        mixdispl = mapping.pure2mix(displcells[d])
        cellfile = chem+'_'+str(d)+'.cell'
        mapping.casprint(mixdispl, cellfile)

def calc_phonons(casfile, supercell='Gamma', method=2, postol=0.001,
                 bands=None, verbose=False):
    """ Compute phonons from completed singlepoint CASTEP calculations of
    perturbations about a relaxed cell.
    Note that this function assumes that the perturbed files are in the
    same directory and use the same naming convention as the parent.

    str casfile : filename (incl. path) to .castep file to compute phonons of
    np.array(3, 3) supercell : size of real space supercell dictates k-points
    If supercell is 'Gamma', only unit cell displacements are created.
    int method : can be 1 (regenerate displacements and assume the same) or
    2 (process all castep files that follow displacement naming convention)
    float postol : Tolerance (Ang) when ion considered displaced (2 only)
    np.array() bands : k-points to compute frequencies at
    bool verbose : function can take a while to run -> prints checkpoints.
    
    returns
    tuple (q_points, distances, frequencies, eigvecs) : output of
    phonopy_instance.get_band_structure() command.
    These are all lists of lists of np.arrays, where the relavent entry of
    frequencies is a list of freqs (in cm-1) of the form np.array(3*phonions)
    and the relevant entry of eigvecs is a matrix of phonon eigenvectors of 
    the form np.array(3*phonions, 3*phonions). """
    # Set up the PHONOPY object for the relaxed cell
    chem = casfile.replace('.castep', '')
    cas = rc.readcas(casfile)
    print("\nreading "+casfile+"\n")
    mixatoms = cas.extract_struc()
    mixkey = cas.get_mixkey()
    mapping = mixmap.mixmap(mixatoms, mixkey)
    pureatoms = mapping.mix2pure(mixatoms)
    pureions = pureatoms.get_number_of_atoms()
    if supercell == 'Gamma':  # Gamma-point means only single unit cell
        supercell = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        bands = [np.array([0.0])]
    PhononObj = phonopy.Phonopy(pureatoms, supercell,
                                factor=phonopy.units.VaspToCm)
    
    if method == 1:
        """ Generate more displs and ASSUME that they are the same as before.
        Less stable than method 2 and should only be used as a check. """
        # Generate perturbed pure cells
        PhononObj.generate_displacements(distance=0.05)
        displcells = PhononObj.get_supercells_with_displacements()
        Ndispls = len(displcells)
        sets_of_forces = np.zeros((Ndispls, pureions, 3))
        for d in range(Ndispls):
            displcasfile = chem+'_'+str(d)+'.castep'
            if verbose:
                print("reading "+displcasfile)
            displmix = rc.readcas(displcasfile)
            displmixforces = displmix.get_forces()
            displpureforces = mapping.mix2pure_forces(displmixforces)
            sets_of_forces[d] = displpureforces
            PhononObj.set_forces(sets_of_forces)
    
    elif method == 2:
        """ Load all the displaced structures and create a displacement_dataset
        object to load to the PHONOPY instance. """
        # Detect displacement files
        allcasfiles = glob.glob(chem+'_*.castep')
        displcasfiles = []
        for casfile in allcasfiles:
            n = casfile.replace(chem+'_', '').replace('.castep', '')
            try:
                int(n)
                displcasfiles.append(casfile)
            except ValueError:
                pass  # file doesn't follow displacement naming convention
        Ndispls = len(displcasfiles)
        first_atoms_list = []
        pureposns = pureatoms.get_positions()
        
        for d in range(Ndispls):
            # Load displacement files and compute displ relative to parent
            displcasfile = chem+'_'+str(d)+'.castep'
            if verbose:
                print("reading "+displcasfile)
            displmix = rc.readcas(displcasfile)
            displmixatoms = displmix.extract_struc()
            displpureatoms = mapping.mix2pure(displmixatoms)
            displpureatoms = pinpos.reorder_atoms(pureatoms, displpureatoms)
            puretotaldispl = displpureatoms.get_positions() - pureposns
            
            # Check that exactly one ion is displaced per perturbed cell
            displions = []
            for i in range(pureions):
                displmagnitude = np.linalg.norm(puretotaldispl[i, :])
                if displmagnitude > postol:
                    displions.append(i)
            if len(displions) > 1:
                raise IndexError('The following ions were all displaced by '
                                 + 'more than ' + str(postol) + ': '
                                 + ' '.join([str(i) for i in displions]))
            elif not displions:
                raise IndexError('No ions were found to be displaced by '
                                 + 'more than ' + str(postol))
            displion = displions[0]
            
            # Extract (pure) forces and construct first_atoms_list
            # this list is later fed to construct displacement_dataset
            displmixforces = displmix.get_forces()
            displpureforces = mapping.mix2pure_forces(displmixforces)
            first_atoms_dict = {}
            first_atoms_dict['number'] = displion
            first_atoms_dict['displacement'] = puretotaldispl[displion, :]
            first_atoms_dict['forces'] = displpureforces
            first_atoms_list += [first_atoms_dict]
            
        # Now construct displacement_dataset (for all perturbations)
        displacement_dataset = {}
        displacement_dataset['natom'] = pureions
        displacement_dataset['first_atoms'] = first_atoms_list
        PhononObj.set_displacement_dataset(displacement_dataset)
    # End of Method 2
    
    # Forces & displacements are now set so may compute the force constants
    if verbose:
        print("\nConstructing force constant matrix")
    PhononObj.produce_force_constants()
    phonopy.harmonic.force_constants.set_translational_invariance(
        PhononObj.force_constants)
    FCs = PhononObj.force_constants
    (FC1, FC2, FC3, FC4) = FCs.shape
    PhononObj.set_force_constants(FCs)
    # Force constants at this point ALREADY obey the ASR near perfectly
    # (i.e. columns/rows sum to zero)
    
    # Compute frequencies from dynamical matrix
    if verbose:
        print("\nDiagonalising dynamical matrix")
    PhononObj.set_band_structure(bands, is_eigenvectors=True)
    return PhononObj.get_band_structure()

