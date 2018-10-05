#!/usr/bin/env python
import argparse
import symeess

parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# Shape input flags
group_symgroup = parser.add_argument_group('Symgroup', 'Symgroup_options')
group_symgroup.add_argument('-m', '--measure', action='store_true', default=False,
                            help='Shape measure of input structure with reference polyhedra')
group_symgroup.add_argument('-c', action='store', type=int, default=None,
                            help='Position of the central atom if exist')
group_symgroup.add_argument('-label', dest='symmetry_operation', action='store',  default=None,
                            help='compute the symmetry operation for the given structure')

# args = parser.parse_args()
args = parser.parse_args(['-m', '-label', 'i',
                          '-o', 'coord',
                          '-input_file', '../examples/ethane.xyz'])

geometry = symeess.file_io.read_input_file(args.input_file)
central_atom = args.c
if args.output_name is not None:
    output_name = args.output_name
else:
    output_name = 'symeess'

example = symeess.Symeess()
example.set_molecules(geometry)
example.write_symgroup_measure(args.symmetry_operation, output_name=output_name)
