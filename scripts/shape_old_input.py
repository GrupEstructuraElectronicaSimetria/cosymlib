#!/usr/bin/env python
import argparse
import symeess
from symeess.shape import shape_tools


parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# Shape input flags
group_shape = parser.add_argument_group('Shape', 'Shape_options')
# group_shape.add_argument('-old_input', action='store_true', default=False)
group_shape.add_argument('-m', '--measure', action='store_true', default=False,
                         help='Shape measure of input structure with reference polyhedra')
group_shape.add_argument('-s', '--structure', action='store_true', default=False,
                         help='Calculate the ideal structure for the input structure')
group_shape.add_argument('-t', '--test', action='store_true', default=False,
                         help='Print the reference structure of the given label')
group_shape.add_argument('-n', action='store', type=str, default=None,
                         help='Print all the possible reference structures of n vertices')

args = parser.parse_args(['-m',
                          '-s',
                          '-input_file', '../examples/coord.dat'])


# Reading and initializing
geometries = symeess.file_io.read_old_input(args.input_file)
geometries, options = geometries
reference_polyhedra = []
for reference in options[1]:
    reference_polyhedra.append(shape_tools.get_shape_label(int(reference), int(options[0][0])))
central_atom = int(options[0][1])
if central_atom == 0:
    central_atom = None

if args.output_name is not None:
    output_name = args.output_name
else:
    output_name = 'symeess'

example = symeess.Symeess()
example.set_molecules(geometries)

# Shape's commands
if args.structure:
    example.write_shape_structure_2file(reference_polyhedra, central_atom=central_atom, output_name=output_name)
if args.measure:
    example.write_shape_measure_2file(reference_polyhedra, central_atom=central_atom, output_name=output_name)
if args.test:
    for reference in reference_polyhedra:
        print(reference)
        print(symeess.shape.shape_tools.get_test_structure(reference, central_atom=central_atom))
if args.n:
    args = parser.parse_args(['-n', '4'])
    print(symeess.shape.shape_tools.get_structure_references(args.n))

# needs parsing
example.write_path_parameters_2file('SP-4', 'T-4', central_atom=central_atom, output_name=output_name)
symeess.write_minimum_distortion_path_shape_2file(reference_polyhedra[0],
                                                  reference_polyhedra[1],
                                                  central_atom=central_atom,
                                                  num_points=20, output_name=output_name)
