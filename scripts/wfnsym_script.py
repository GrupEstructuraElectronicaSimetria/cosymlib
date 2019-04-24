#!/usr/bin/env python
import argparse
from symeess import Symeess, file_io


parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# args = parser.parse_args()
args = parser.parse_args(['-o', 'TiH4',
                          '-input_file', '/Users/efrem/PycharmProjects/ExtendedHuckel/test_Tih4.fchk'])
                          #'-input_file', '../examples/pirrol.in.fchk'])


molecule = file_io.get_molecule_from_file_fchk(args.input_file)
if args.output_name is not None:
    output_name = args.output_name
else:
    output_name = 'symeess'

example = Symeess()
example.set_molecules(molecule)
axis1 = [1.000000000, 1.000000000, 1.000000000]
axis2 = [1.000000000, -1.000000000, -1.000000000]
center_operation = [0., 0., 0.]
example.write_wnfsym_measure_2file('Td', vector_axis1=axis1,
                                   vector_axis2=axis2,
                                   center=center_operation,
                                   output_name=output_name)
