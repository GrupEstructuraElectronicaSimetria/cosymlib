#!/usr/bin/env python
import argparse
from symeess import Symeess, file_io

parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# args = parser.parse_args()
args = parser.parse_args(['-input_file' , '../examples/CpTiCl3_cart.fchk'])

molecule = file_io.read_input_file(args.input_file)
print(molecule)
quit()
symeess = Symeess(molecule)
symeess.write_wnfsym_measure_2file('Oh', [-2.027247, 0.000133, -0.898469], [0.40757934076903307, 1.746331, -0.919377],
                                   [0.002440,   -0.000122, 0.017307])
