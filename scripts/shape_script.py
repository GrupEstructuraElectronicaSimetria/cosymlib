#!/usr/bin/env python
import argparse
import symeess

parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# Shape input flags
group_shape = parser.add_argument_group('Shape', 'Shape_options')
group_shape.add_argument('--old_input', action='store_true', default=False)
group_shape.add_argument('-m', '--measure', action='store_true', default=False,
                         help='Shape measure of input structure with reference polyhedra')
group_shape.add_argument('-s', '--structure', action='store_true', default=False,
                         help='Calculate the ideal structure for the input structure')
group_shape.add_argument('-t', '--test', action='store_true', default=False,
                         help='Print the reference structure of the given label')
group_shape.add_argument('-c', action='store', type=int, default=None,
                         help='Position of the central atom if exist')
group_shape.add_argument('-label', dest='reference_polyhedron', action='store',  default=None,
                         help='Use labels from Shape manual for desire reference polyhedra')
group_shape.add_argument('-n', action='store', type=str, default=None,
                         help='Print all the possible reference structures of n vertices')

# args = parser.parse_args(['-n', '4'])
# print(shape.get_structure_references(args.n))
args = parser.parse_args(['-m', '-label', 'SP-4 T-4',
                          '-c', '1',
                          '-o', '../examples/coord',
                          '-input_file', '../examples/coord.xyz'])


reference_polyhedron = args.reference_polyhedron.split()
geometries = symeess.file_io.read_input_file(args.input_file)

example = symeess.Symeess()
example.set_molecules(geometries)

# Shape's commands
if args.structure:
    example.write_shape_structure_2file(reference_polyhedron, central_atom=central_atom)
if args.measure:
    example.write_shape_measure_2file(reference_polyhedron, central_atom=central_atom)
if args.test:
    for reference in reference_polyhedron:
        print(reference)
        print(shape_tools.get_test_structure(reference, central_atom=central_atom))
if args.n:
    args = parser.parse_args(['-n', '4'])
    print(shape_tools.get_structure_references(args.n))

example.write_path_parameters_2file('SP-4', 'T-4' , central_atom=central_atom)
symeess.write_minimum_distortion_path_shape_2file(reference_polyhedron[0],
                                                  reference_polyhedron[1],
                                                  central_atom=central_atom,
                                                  num_points=20)

example = symeess.Symeess()
example.set_molecules(geometries)
example.write_shape_structure_2file(reference_polyhedron, central_atom=args.c)
example.write_shape_measure_2file(reference_polyhedron, central_atom=args.c)
example.write_path_parameters_2file('SP-4', 'T-4' , central_atom=args.c)
symeess.write_minimum_distortion_path_shape_2file(reference_polyhedron[0],
                                                  reference_polyhedron[1], central_atom=args.c)
