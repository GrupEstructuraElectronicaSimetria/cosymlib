#!/usr/bin/env python
import argparse
from symeess import Symeess, file_io
import numpy as np


parser = argparse.ArgumentParser(description='Symeess')

parser.add_argument('-input_file', type=str, help='input file name(+extension)')
parser.add_argument('-o', '-output', dest='output_name', default=None, help='output_name')

# args = parser.parse_args()
args = parser.parse_args(['-input_file' , '../examples/CpTiCl3.fchk'])

molecule = file_io.read_input_file(args.input_file)
example = Symeess()
example.set_molecules(molecule)
VAxis1 = [0.000000000, 0.000000000, 1.000000000]
VAxis2 = [-2.027247, 0.000133, -0.898469]
# norm = np.linalg.norm(VAxis2)
# VAxis2 = VAxis2/norm
RCread = [0.002440, -0.000122, 0.017307]
example.write_wnfsym_measure_2file('Td', vector_axis1=VAxis1,
                                   vector_axis2=VAxis2,
                                   center=RCread)
