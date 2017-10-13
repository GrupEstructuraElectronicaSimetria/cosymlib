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
        name = lines.readline()
        input_molecule = lines.read().splitlines()
        molecule = Molecule(structure=input_molecule)

        return molecule


def read_file_cor(file_name):
    with open(file_name, mode='r') as lines:
        input_molecule = []
        molecules = {}

        name = lines.readline().split()[0]
        for line in lines:
            try:
                float(line.split()[1])
                input_molecule.append(line.split()[:-1])
            except ValueError:
                molecules[name] = Molecule(structure=input_molecule)
                input_molecule = []
                name = line.split()[0]
        molecules[name] = Molecule(structure=input_molecule)

        return molecules


def write(output_name, data, shape_label=None, shape_choices=None):

    if shape_choices:
        output = open(output_name, 'w') if output_name else sys.stdout
        output.write('-'*80 + '\n')
        output.write('Shape measure/s with {} reference \n'.format(shape_label))
        output.write('-'*80 + '\n')
        for element in shape_choices:
            output.write('{}\n'.format(' '.join(element.split('_')[1:])))
            for key in data:
                if isinstance(data[key][element], np.ndarray):
                    output.write('{} \n'.format(key))
                    for item in data[key][element]:
                        output.write('{:11.7f} {:11.7f} {:11.7f}\n'.format(item[0], item[1], item[2]))
                else:
                    output.write('{} {:4.3f}\n'.format(key, float(data[key][element])))
            output.write('\n')
