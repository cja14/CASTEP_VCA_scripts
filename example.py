import readmixcastep as rc
import mixmap

""" A silly script to test that all the functionality works """

########################################################
"""
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
"""
########################################################

# Read pure cell

print('\n\n\nPure Cell\n\n\n')

cas = rc.readcell('examples/Ca2.15Sr0.85Ti2O7_Amam_nonSS.cell')

mixkey = {'Ca': {'Ca': 0.7166, 'Sr': 0.2834},
          'Ti': {'Ti': 1.0},
          'O': {'O': 1.0}}

atoms = cas.extract_struc()

mixatoms = mixmap.mixmap.mix_atoms(atoms, mixkey)
mapping = mixmap.mixmap(mixatoms, mixkey)

mapping.casprint(mixatoms, 'test3.cell', pure=False)  # write mix cell
