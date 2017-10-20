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


# def write(output_name, data):
#
#     output = open(output_name, 'w') if output_name else sys.stdout
#
#     output.write('-'*40 + '\n')
#     output.write('{}'.format('END of CALCULATION'.rjust(26)))
#     output.close()


def write_shape_data(output_name, data, shape_label, molecule_names, option):

    extensions = {'measure': '.tab', 'structure': '.out', 'test': '.tst'}
    output = open('../examples/'+output_name + extensions[option], 'w')

    output.write('-'*40 + '\n')
    output.write('Shape measure/s \n')
    output.write('-'*40 + '\n')

    if 'measure' in option:
        output.write("{} {}\n".format('measure'.upper(), '    '.join(shape_label).rjust(14)))
        output.write('\n')
        for idx, molecule in enumerate(data):
            output.write('{}'.format(molecule_names[idx]))
            if molecule_names[idx].strip() == '':
                n = 4 + len(molecule_names[idx])
            else:
                n = 14 - len(molecule_names[idx])
            for label in shape_label:
                output.write(' {:{width}.{prec}f}'.format(molecule[label]['measure'], width=n, prec=3))
                n = 7
            output.write('\n')
        output.write('\n')

    if 'structure' in option:
        output.write("{}\n".format('ideal_structure'.upper()))
        for idx, molecule in enumerate(data):
            output.write('\n')
            output.write('{}'.format(molecule_names[idx]))

            n = int(23 - len(molecule_names[idx]))
            for label in shape_label:
                output.write('{}'.format(label.rjust(n)))
                n = 36 + len(label)
            output.write('\n')

            for idn, symbol in enumerate(molecule['symbols']):
                output.write('{:3s}'.format(molecule['symbols'][idn]))
                for label in shape_label:
                    array = molecule[label]['structure'][idn]
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

    # output.write('-'*40 + '\n')
    # output.write('-'*40 + '\n')
    output.close()
