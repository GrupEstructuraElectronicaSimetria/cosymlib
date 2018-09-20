import symeess.file_io as file_io
import symeess.shape as shape


def print_data(geometries):
    print('{:10} {:^10} {:^10} {:^10}'.format('name', 'SP-4', 'SS-4', 'PP-5'))
    print('-'*43)
    for geometry in geometries:
        print('{:10} {:^10.3f} {:^10.3f} {:^10.3f}'.format(geometry.get_name(),
                                                           geometry.get_shape_measure('SP-4', central_atom=1),
                                                           geometry.get_shape_measure('SS-4', central_atom=1),
                                                           geometry.get_shape_measure('PP-5')
                                                           ))
    print()


# Get structures from files
molecules_set = file_io.read_geometry_from_xyz_file('coord.xyz')
fragments_set = file_io.read_geometry_from_cor_file('coord.cor')

# Call shape as method of Geometry class
print_data(molecules_set)
print_data(fragments_set)

# Check multiple calls of shape one calculatioin
methane = molecules_set[0]
for i in range(100):
    measure = methane.get_shape_measure('SP-4', central_atom=1)

print('final measure:', measure)

# Call shape as method of Shape class (semi function call)
print('measure:', shape.Shape(methane).measure('SP-4', central_atom=1))
print('structure:\n', shape.Shape(methane).structure('SP-4', central_atom=1))
