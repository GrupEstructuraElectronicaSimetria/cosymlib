import os
import sys
import numpy as np
from molecule import Molecule


def read(input_name):

    file_name, file_extension = os.path.splitext(input_name)
    method_name = 'read_file_' + file_extension[1:]
    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(method_name)
    if not method:
        raise NotImplementedError("Method %s not implemented" % method_name)

    return method(input_name)


def read_file_xyz(file_name):
    with open(file_name, mode='r') as lines:
        lines.readline()
        lines.readline()
        input_molecule = lines.read().splitlines()
        molecule = Molecule(structure=input_molecule)

        return molecule


def read_file_cor(file_name):
    with open(file_name, mode='r') as lines:
        input_molecule = []
        lines.readline()
        for line in lines:
            input_molecule.append(line.split()[:-1])
        molecule = Molecule(structure=input_molecule)

        return molecule


def write(output_name, data, shape_label=None, shape_choices=None):

    if shape_choices:
        output = open(output_name, 'w') if output_name else sys.stdout
        output.write('-'*80 + '\n')
        output.write('Shape measures with {} reference \n'.format(shape_label))
        output.write('-'*80 + '\n')
        for element in shape_choices:
            output.write('{}\n'.format(' '.join(element.split('_')[1:])))
            if isinstance(data[element], np.ndarray):
                for item in data[element]:
                    output.write('{:11.7f} {:11.7f} {:11.7f}\n'.format(item[0], item[1], item[2]))
            else:
                output.write('{:3.2f}\n'.format(float(data[element])))
