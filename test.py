import readmixcastep as rc
import mixmap
import mixcastep as mc

cas = rc.readcell('examples/Ca2.15Sr0.85Ti2O7_Amam.cell')
mixatoms = cas.extract_struc()
mixkey = cas.get_mixkey()

casatoms = mc.casatoms(mixatoms, mixkey=mixkey)

print(casatoms)

# Calculation params

casatoms.read_cellparams('examples/Ca2.15Sr0.85Ti2O7_Amam.cell')

print('\nkpoints_mp_grid:')
print(casatoms.get_keyword('kpoints_mp_grid'))

casatoms.set_keyword('kpoints_mp_grid', '2  2  1')

print('\nsymmetry_generate:')
print(casatoms.get_keyword('symmetry_generate'))

casatoms.set_keyword('symmetry_generate', False)

print('\ncell_constraints:')
print(casatoms.get_block('cell_constraints'))

casatoms.set_block('cell_constraints', [[1, 2, 3], [4, 5, 6]])

print('\nspecies_pot:')
print(casatoms.get_block('species_pot'))

casatoms.set_block('species_pot', False)

# Write cell file
casatoms.writecell('test.cell')

"""
mapping = mixmap.mixmap(mixatoms, mixkey)
pureatoms = mapping.mix2pure(mixatoms)
mixatoms2 = mixmap.create_mixture(pureatoms, mixkey)
elem_mixkey = mapping.site2elem_mixkey(mixkey)
mixatoms3 = mixmap.create_mixture(pureatoms, elem_mixkey)


print('\nmixatoms:\n'+str(mixatoms))
print('\npureatoms:\n'+str(pureatoms))
print('\nmixatoms2:\n'+str(mixatoms2))
print('\nmixatoms3:\n'+str(mixatoms3))
"""
