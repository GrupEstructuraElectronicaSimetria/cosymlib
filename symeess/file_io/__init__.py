import os
import sys
from itertools import islice
from symeess.molecule import Molecule


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
        molecules = []

        while True:
            try:
                input_molecule = []
                n_atoms = int(lines.readline())
                name = lines.readline()
                if name.strip() != '':
                    name = name.split()[0]
                for line in list(islice(lines, n_atoms)):
                    input_molecule.append(line.split())
                molecules.append(Molecule(input_molecule, name=name))
            except ValueError:
                return molecules


def read_file_cor(file_name):
    with open(file_name, mode='r') as lines:
        input_molecule = []
        molecules = []
        name = lines.readline().split()[0]
        for line in lines:
            try:
                float(line.split()[1])
                input_molecule.append(line.split()[:-1])
            except ValueError:
                molecules.append(Molecule(input_molecule, name=name))
                input_molecule = []
                name = line.split()[0]
        molecules.append(Molecule(input_molecule, name=name))

        return molecules


def write_shape_data(data, shape_label, molecule_names, option, output_name=sys.stdout):

    extensions = {'measure': '.tab', 'structure': '.out', 'test': '.tst'}
    if not os.path.exists('./results'):
        os.makedirs('./results')
    output = open('results/'+output_name + extensions[option], 'w')

    output.write('-'*40 + '\n')
    output.write('Shape measure/s \n')
    output.write('-'*40 + '\n')

    if 'measure' in option:
        output.write('{}'.format('measure'.upper()))
        for label in shape_label:
            n = len(label) + 4
            output.write('{}'.format(label.rjust(n)))
        output.write('\n')
        for idx, molecule_name in enumerate(molecule_names):
            output.write('{}'.format(molecule_name))
            if molecule_names[idx].strip() == '':
                n = 4 + len(molecule_names[idx])
            else:
                n = 14 - len(molecule_names[idx])
            for label in shape_label:
                output.write(' {:{width}.{prec}f}'.format(data[label][idx], width=n, prec=3))
                n = 7
            output.write('\n')
        output.write('\n')

    if 'structure' in option:
        output.write("{}\n".format('ideal_structure'.upper()))
        for idx, molecule_name in enumerate(molecule_names):
            output.write('\n')
            output.write('{}'.format(molecule_names[idx]))

            n = int(23 - len(molecule_names[idx]))
            for label in shape_label:
                output.write('{}'.format(label.rjust(n)))
                n = 37
            output.write('\n')

            for idn, symbol in enumerate(data['symbols'][idx]):
                output.write('{:2s}'.format(symbol))
                for label in shape_label:
                    array = data[label][idx][idn]
                    output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
                output.write('\n')

        output.write('\n')

    if 'test_structure' in option:
        output.write("{}\n".format('test_structure'.upper()))
        n = 20
        for label in shape_label:
            output.write('{}'.format(label.rjust(n)))
            n = 36 + len(label)
        output.write('\n')

        for idx in list(range(len(data[0]['symbols']))):
            for label in shape_label:
                array = data[0][label]['test_structure'][idx]
                output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
            output.write('\n')

    output.close()


def write_shape_map_2file(shape_label1, shape_label2, path, output_name='symeess_shape_map'):
    output = open('results/' + output_name + '.pth', 'w')
    # output.write("{}\n".format("minimal distortion path"))
    output.write(" {:6} {:6}\n".format(shape_label1, shape_label2))
    for idx, value in enumerate(path[0]):
        output.write('{:6.3f} {:6.3f}'.format(path[0][idx], path[1][idx]))
        output.write('\n')
