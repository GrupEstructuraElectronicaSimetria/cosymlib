#!/usr/bin/env python
import argparse
from symeess import Symeess, file_io, shape

parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
# parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# Shape input flags
group_shape = parser.add_argument_group('Shape', 'Shape_options')
group_shape.add_argument('-old_input', action='store_true', default=False)
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
# print(shape.get_shape_references(args.n))
args = parser.parse_args(['-m', '-label', 'SP-4 T-4',
                          '-c', '1',
                          # '-o', '../examples/coord.tab',
                          '-old_input',
                          '-input_file' , '../examples/coord.dat'])

molecules = file_io.read(args.input_file, args.old_input)
if args.old_input:
    molecules, options = molecules
    central_atom = options[0][1]
    if central_atom == 0:
        central_atom = None
    reference_polyhedra = ''
    for reference in options[1]:
        reference_polyhedra += shape.get_shape_label(int(reference), int(options[0][0]))+' '
else:
    reference_polyhedra = args.reference_polyhedra
symeess = Symeess(molecules)
symeess.write_shape_structure_2file(reference_polyhedra, central_atom=args.c)
symeess.write_shape_measure_2file(reference_polyhedra, central_atom=args.c)
symeess.write_path_parameters_2file('SP-4', 'T-4' , central_atom=args.c)
symeess.write_minimum_distortion_path_shape_2file('SP-4', 'T-4', central_atom=args.c, num_points=50)