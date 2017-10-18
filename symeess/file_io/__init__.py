import os
import sys
from itertools import islice
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
        molecules = {}
        names = []

        while True:
            try:
                input_molecule = []
                n_atoms = int(lines.readline())
                name = lines.readline().split()[0]
                names.append(name)
                for line in list(islice(lines, n_atoms)):
                    input_molecule.append(line.split())
                molecules[name] = Molecule(structure=input_molecule)
            except ValueError:
                return molecules, names


def read_file_cor(file_name):
    with open(file_name, mode='r') as lines:
        input_molecule = []
        molecules = {}
        names = []
        name = lines.readline().split()[0]
        names.append(name)
        for line in lines:
            try:
                float(line.split()[1])
                input_molecule.append(line.split()[:-1])
            except ValueError:
                molecules[name] = Molecule(structure=input_molecule)
                input_molecule = []
                name = line.split()[0]
                names.append(name)
        molecules[name] = Molecule(structure=input_molecule)

        return molecules, names


# def write(output_name, data):
#
#     output = open(output_name, 'w') if output_name else sys.stdout
#
#     output.write('-'*40 + '\n')
#     output.write('{}'.format('END of CALCULATION'.rjust(26)))
#     output.close()


def write_shape(output_name, data, shape_label, molecule_names):

    output = open(output_name, 'w') if output_name else sys.stdout
    options = sorted(list(data[molecule_names[0]][shape_label[0]].keys()))

    output.write('-'*40 + '\n')
    output.write('Shape measure/s \n')
    output.write('-'*40 + '\n')

    if 'measure' in options:
        output.write("{} {}\n".format('measure'.upper(),'   '.join(shape_label)))
        for molecule in molecule_names:
            output.write('{}'.format(molecule))
            for label in shape_label:
                output.write(' {:4.3f}'.format(data[molecule][label]['measure']))
            output.write('\n')
        output.write('\n')

    if 'ideal_structure' in options:
        output.write("{}\n".format('ideal_structure'.upper()))
        for molecule in molecule_names:
            output.write('\n')
            output.write('{}'.format(molecule))
            n = 0
            for label in shape_label:
                n += 15
                output.write('{}'.format(label.rjust(n)))
            output.write('\n')

            for idx, symbol in enumerate(data[molecule]['symbols']):
                output.write('{:3s}'.format(symbol))
                for label in shape_label:
                    array = data[molecule][label]['ideal_structure'][idx]
                    output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
                output.write('\n')

        output.write('\n')

    if 'test_structure' in options:
        output.write("{}\n".format('test_structure'.upper()))
        for molecule in molecule_names:
            output.write('\n')
            output.write('{}'.format(molecule))
            n = 0
            for label in shape_label:
                n += 15
                output.write('{}'.format(label.rjust(n)))
            output.write('\n')

            for idx, symbol in enumerate(data[molecule]['symbols']):
                output.write('{:3s}'.format(symbol))
                for label in shape_label:
                    array = data[molecule][label]['test_structure'][idx]
                    output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
                output.write('\n')

    output.write('-'*40 + '\n')
    output.write('-'*40 + '\n')
    output.close()
