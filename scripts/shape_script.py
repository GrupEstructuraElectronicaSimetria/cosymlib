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
group_shape.add_argument('-label', dest='reference_polyhedra', action='store',  default=None,
                         help='Use labels from Shape manual for desire reference polyhedra')
group_shape.add_argument('-n', action='store', type=str, default=None,
                         help='Print all the possible reference structures of n vertices')

# args = parser.parse_args(['-n', '4'])
# print(shape.get_structure_references(args.n))
args = parser.parse_args(['-m', '-label', 'SP-4 T-4',
                          '-c', '1',
                          '-o', '../examples/coord',
                          '-input_file' , '../examples/coord.xyz'])

molecules = symeess.file_io.read_input_file(args.input_file)
reference_polyhedra = args.reference_polyhedra.split()
example = symeess.Symeess()
example.set_molecules(molecules)
example.write_shape_structure_2file(reference_polyhedra, central_atom=args.c)
example.write_shape_measure_2file(reference_polyhedra, central_atom=args.c)
example.write_path_parameters_2file('SP-4', 'T-4' , central_atom=args.c)
symeess.write_minimum_distortion_path_shape_2file(reference_polyhedra[0],
                                                  reference_polyhedra[1], central_atom=args.c)
