import readmixcastep as rc
import mixmap
import phonons_VCA as pVCA
import os

""" A silly script to test that all the functionality works.
Running this script should test that most functions and classes run without
errors (although as yet no checks test that the results are sensible). """

########################################################
########################################################

###################################
# Manipulating pure & mixed cells #
###################################

# Read mix cell

print('\n\n\nMix Cell\n\n\n')

cas = rc.readcell('examples/Ca2.15Sr0.85Ti2O7_Amam.cell')

print(cas.get_kpoints())
print(cas.get_psps())
print(cas.get_elements())
print(cas.get_init_spin())
print(cas.get_mixkey())
print(cas.get_posns())
print(cas.get_ext_press())
print(cas.get_cell_constrs())
print(cas.get_cell())

mixkey = cas.get_mixkey()

########################################################

# Read mix castep

print('\n\n\nMix Castep\n\n\n')

cas = rc.readcas('examples/Ca2.15Sr0.85Ti2O7_Amam.castep')

print(cas.get_kpoints())
print(cas.get_psps())
print(cas.get_elements())
print(cas.get_init_spin())
print(cas.get_mixkey())
print(cas.get_posns())
print(cas.get_ext_press())
print(cas.get_cell_constrs())
print(cas.get_cell())
print(cas.get_enthalpy())
print(cas.get_energy())
print(cas.get_forces())
print(cas.get_stresses())
print(cas.get_Fmax())

mixkey = cas.get_mixkey()

atoms = cas.extract_struc()
mapping = mixmap.mixmap(atoms, mixkey)

print(mapping.pure2mix_map)
print(mapping.mix2pure_map)
print(mapping.mixsitemixes)

print(mapping.mixions)

mapping.casprint(atoms, 'test1.cell', pure=False)  # write mix cell
mapping.casprint(atoms, 'test2.cell', pure=True)  # write pure cell

########################################################

# Test SS_to_endmember.py command line tool

# On a .castep file

os.system("python SS_to_endmember.py examples/Ca2.15Sr0.85Ti2O7_Amam.castep")

print("\nRan SS_to_endmember.py script on .castep file.")

# And on a .cell file

os.system("python SS_to_endmember.py examples/Ca2.15Sr0.85Ti2O7_Amam.cell")

print("\nRan SS_to_endmember.py script on .cell file.")

########################################################

# Read pure cell

print('\n\n\nPure Cell\n\n\n')

cas = rc.readcell('examples/Ca2.15Sr0.85Ti2O7_Amam_nonSS.cell')

mixkey = {'Ca': {'Ca': 0.7166, 'Sr': 0.2834},
          'Ti': {'Ti': 1.0},
          'O': {'O': 1.0}}

atoms = cas.extract_struc()

mixatoms = mixmap.create_mixture(atoms, mixkey)
mapping = mixmap.mixmap(mixatoms, mixkey)

mapping.casprint(mixatoms, 'test3.cell', pure=False)  # write mix cell

########################################################
########################################################

######################
# Phonon calcuations #
######################

# Generate perturbations for phonon calcuation

pVCA.gen_perturbations('examples/Ca1.5Sr0.5GeO4/Ca1.5Sr0.5GeO4.castep')

print('\nGenerated perturbed cell files for phonon input.')

########################################################

# Compute phonons (using method 1)

output = pVCA.calc_phonons('examples/Ca1.5Sr0.5GeO4/Ca1.5Sr0.5GeO4.castep',
                           method=1, verbose=True)
q_points, distances, frequencies, eigenvectors = output
freqs = frequencies[0][0]
eigvecs = eigenvectors[0][0]

print("\nFirst 10 frequencies:")
print(freqs[:10])

print("\nFirst 5 displacments of first two eigenvectors:")
print(eigvecs[1, :15].reshape(5, 3))
print(eigvecs[2, :15].reshape(5, 3))

########################################################

# Compute phonons (using method 2) -- slower but more stable method

output = pVCA.calc_phonons('examples/Ca1.5Sr0.5GeO4/Ca1.5Sr0.5GeO4.castep',
                           method=2, verbose=True)
q_points, distances, frequencies, eigenvectors = output
freqs = frequencies[0][0]
eigvecs = eigenvectors[0][0]

print("\nFirst 10 frequencies:")
print(freqs[:10])

print("\nFirst 5 displacments of first two eigenvectors:")
print(eigvecs[1, :15].reshape(5, 3))
print(eigvecs[2, :15].reshape(5, 3))

########################################################

print("\n\nAll functions and methods appeared to run succesfully.\n\n")
