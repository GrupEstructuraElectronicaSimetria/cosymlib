import cosymlib.file_io as file_io
import cosymlib.shape as shape
from cosymlib.symmetry import Symmetry
from cosymlib import Cosymlib


def print_shape_data(geometries):
    print('{:10} {:^10} {:^10} {:^10}'.format('name', 'SP-4', 'SS-4', 'PP-5'))
    print('-'*43)
    for geometry in geometries:
        print('{:10} {:^10.3f} {:^10.3f} {:^10.3f}'.format(geometry.name,
                                                           geometry.get_shape_measure('SP-4', central_atom=1),
                                                           geometry.get_shape_measure('SS-4', central_atom=1),
                                                           geometry.get_shape_measure('PP-5')
                                                           ))
    print()


def print_csm(data):
    print('\nWaveFunction: CSM-like values')
    print('     ' + '  '.join(['{:^7}'.format(s) for s in data['labels']]))
    print('Grim' + '  '.join(['{:7.3f}'.format(s) for s in data['grim']]))
    print('CSM ' + '  '.join(['{:7.3f}'.format(s) for s in data['csm']]))


# Get structures from files
geometries_list = file_io.get_geometry_from_file_xyz('coord.xyz', read_multiple=True)
fragments_list = file_io.get_geometry_from_file_cor('coord.cor', read_multiple=True)

# Call shape as method of Geometry class
print_shape_data(geometries_list)
print_shape_data(fragments_list)


# Check multiple calls of shape one calculatioin
methane = geometries_list[0]
for i in range(100):
    measure = methane.get_shape_measure('SP-4', central_atom=1)

print('final measure:', measure)

# Call shape as method of Shape class (semi function call)
print('measure:', shape.Shape(methane).measure('SP-4', central_atom=1))
print('structure:\n', shape.Shape(methane).structure('SP-4', central_atom=1))

# test symgroup
print('\nSYMMETRY\n--------')

geometries_list = file_io.get_geometry_from_file_pdb('methane.pdb', read_multiple=True)
print('measure C3: {} '.format(geometries_list[0].get_symmetry_measure('c3')))

# Check multiple calls of symgroup one calculatioin
for i in range(100):
    measure = methane.get_symmetry_measure('C3', central_atom=1)

print('measure: {:^10.3f} '.format(measure))
print('measure: {:^10.3f} '.format(methane.get_symmetry_measure('C4', central_atom=1, fix_permutation=False)))

# Call symgroup as method of Symgroup class (semi function call)
print('measure: {:^10.3f} '.format(Symmetry(methane, central_atom=1).measure('C4')))
print(Symmetry(methane, central_atom=1).nearest_structure('C4'))


# test WFNSYM
print('\nWFNSYM\n--------')


molecule = file_io.get_molecule_from_file_fchk('pirrol.fchk', read_multiple=False)

measure_dict = molecule.get_wf_symmetry('C2v',
                                        axis=[0.000000, 0.000000, 1.000000],  # valor defecte
                                        # vector_axis2=[-2.027247,  0.000133, -0.898469],
                                        center=[0.002440, -0.000122,  0.017307])
print_csm(measure_dict)

print('\nCOSYMLIB\n--------')
geometries_list = file_io.get_geometry_from_file_xyz('coord_2.xyz', read_multiple=True)
molecules_set = Cosymlib(geometries_list)
print('\n\n***********{}************'.format('print_electronic_symmetry_measure()'))
molecules_set.print_electronic_symmetry_measure('C2v', axis=[0, 0, 1], center=[0.0, 0.0, 0.0])


molecules_list = file_io.read_generic_structure_file('pirrol.fchk', read_multiple=True)
molecules_set = Cosymlib(molecules_list)
print('\n\n***********{}************'.format('print_info()'))
molecules_set.print_info()
print('\n\n***********{}************'.format('print_wnfsym_irreducible_repr()'))
molecules_set.print_wnfsym_irreducible_repr('C2v', axis=[0, 0, 1], center=[0.0, 0.0, 0.0])
print('\n\n***********{}************'.format('print_wnfsym_sym_ovelap()'))
molecules_set.print_wnfsym_sym_ovelap('C2v', axis=[0, 0, 1], center=[0.0, 0.0, 0.0])
print('\n\n***********{}************'.format('print_wnfsym_sym_matrices()'))
molecules_set.print_wnfsym_sym_matrices('C2v', axis=[0, 0, 1], center=[0.0, 0.0, 0.0])
