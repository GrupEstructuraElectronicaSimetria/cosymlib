import cosymlib.file_io as file_io
import cosymlib.shape as shape
from cosymlib.symmetry import Symmetry


def print_shape_data(geometries):
    print('{:10} {:^10} {:^10} {:^10}'.format('name', 'SP-4', 'SS-4', 'PP-5'))
    print('-'*43)
    for geometry in geometries:
        print('{:10} {:^10.3f} {:^10.3f} {:^10.3f}'.format(geometry.get_name(),
                                                           geometry.get_shape_measure('SP-4', central_atom=1),
                                                           geometry.get_shape_measure('SS-4', central_atom=1),
                                                           geometry.get_shape_measure('PP-5')
                                                           ))
    print()


def print_csm(data):
    print('\nWaveFunction: CSM-like values')
    print('     ' + '  '.join(['{:^7}'.format(s) for s in data.SymLab]))
    print('Grim' + '  '.join(['{:7.3f}'.format(s) for s in data.grim_coef]))
    print('CSM ' + '  '.join(['{:7.3f}'.format(s) for s in data.csm_coef]))


# Get structures from files
molecules_set = file_io.get_geometry_from_file_xyz('coord.xyz', read_multiple=True)
fragments_set = file_io.get_geometry_from_file_cor('coord.cor', read_multiple=True)

# Call shape as method of Geometry class
print_shape_data(molecules_set)
print_shape_data(fragments_set)


# Check multiple calls of shape one calculatioin
methane = molecules_set[0]
for i in range(100):
    measure = methane.get_shape_measure('SP-4', central_atom=1)

print('final measure:', measure)

# Call shape as method of Shape class (semi function call)
print('measure:', shape.Shape(methane).measure('SP-4', central_atom=1))
print('structure:\n', shape.Shape(methane).structure('SP-4', central_atom=1))

# test symgroup
print('\nSYMMETRY\n--------')

molecule = file_io.get_geometry_from_file_pdb('methane.pdb', read_multiple=False)
print('measure C3: {} '.format(molecule.get_symmetry_measure('c3', center=[0, 0, 0])))

# Check multiple calls of symgroup one calculatioin
for i in range(100):
    measure = methane.get_symmetry_measure('C3', central_atom=1)

print('measure: {:^10.3f} '.format(measure))
print('measure: {:^10.3f} '.format(methane.get_symmetry_measure('C4', central_atom=1)))

# Call symgroup as method of Symgroup class (semi function call)
print('measure: {:^10.3f} '.format(Symmetry(methane, central_atom=1).measure('C4')))
print(Symmetry(methane, central_atom=1).nearest_structure('C4'))


# test WFNSYM
print('\nWFNSYM\n--------')


molecule = file_io.get_molecule_from_file_fchk('pirrol.fchk')
molecules_set = file_io.read_input_file('pirrol.fchk', read_multiple=True)

data = molecule.get_mo_symmetry('C2v',
                                 vector_axis1=[ 0.000000,  0.000000,  1.000000],  # valor defecte
                                 # vector_axis2=[-2.027247,  0.000133, -0.898469],
                                 center=[0.002440, -0.000122,  0.017307])  # valor per defecte (CM)

print_csm(data)
