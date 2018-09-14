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

args = parser.parse_args(['-input_file' , '../examples/CpTiCl3_cart.fchk'])

# args = parser.parse_args()
molecules = file_io.read(args.input_file, args.old_input)
symeess = Symeess(molecules)
symeess.write_wnfsym_measure_2file('Td', [-2.027247, 0.000133, -0.898469], [0.40757934076903307, 1.746331, -0.919377],
                                   [0.002440, -0.000122, 0.017307])
