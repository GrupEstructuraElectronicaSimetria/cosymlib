import os
import sys
from symeess.molecule import Molecule


# INPUT part
def read(input_name, old_input=False):
    if old_input:
        return read_old_input(input_name)
    else:
        file_name, file_extension = os.path.splitext(input_name)
        method_name = 'read_file_' + file_extension[1:]
        possibles = globals().copy()
        possibles.update(locals())
        method = possibles.get(method_name)
        if not method:
            raise NotImplementedError("Method %s not implemented yet" % method_name)
        return method(input_name)


def read_file_fchk(file_name):
    key_list = ['Charge', 'Multiplicity', 'Number of electrons', 'Atomic numbers', 'Current cartesian coordinates',
                'Shell type', 'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha MO coefficients', 'Beta MO coefficients']
    input_molecule = [[] for _ in range(len(key_list))]
    read = False
    with open(file_name, mode='r') as lines:
        name = lines.readline()
        line = lines.readline().split()
        if 'R' in line[1]:
            del key_list[-1]

        n = 0
        options = True
        for line in lines:
            if read:
                try:
                    float(line.split()[0])
                    input_molecule[n].append(line.split())

                except ValueError:
                    input_molecule[n] = reformat_input(input_molecule[n])
                    read = False


            for idn, key in enumerate(key_list):
                if key in line:
                    if n == len(key_list) - 1:
                        break
                    if options and idn != 3:
                        input_molecule[idn].append(int(line.split()[-1]))
                        n = idn
                        break
                    else:
                        options = False
                    if n == idn:
                        n += 1
                    else:
                        n = idn
                    read = True
                    break
        return [Molecule(structure_data=input_molecule[3:5],
                         electronic_structure=input_molecule[:3]+input_molecule[5:],
                         name=name)]


def read_file_xyz(file_name):
    input_molecule = [[], []]
    molecules = []
    with open(file_name, mode='r') as lines:
        name = lines.readline().split()[0]
        lines.readline()
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    float(line.split()[1])
                    input_molecule[0].append(line.split()[0])
                    input_molecule[1].append(line.split()[1:])
                except (ValueError, IndexError):
                    if input_molecule:
                        molecules.append(Molecule(input_molecule, name=name))
                    input_molecule = [[], []]
                    lines.readline()
                    name = line.split()[0]
        molecules.append(Molecule(input_molecule, name=name))
    return molecules


def read_file_cor(file_name):
    input_molecule = [[], []]
    molecules = []
    with open(file_name, mode='r') as lines:
        name = lines.readline().split()[0]
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    float(line.split()[1])
                    input_molecule[0].append(line.split()[0])
                    input_molecule[1].append(line.split()[1:-1])
                except (ValueError, IndexError):
                    if input_molecule:
                        molecules.append(Molecule(input_molecule, name=name))
                    input_molecule = [[], []]
                    name = line.split()[0]
        molecules.append(Molecule(input_molecule, name=name))
    return molecules


def read_old_input(file_name):
    input_molecule = [[], []]
    molecules = []
    options = []
    with open(file_name, mode='r') as lines:
        while len(options) < 2:
            line = lines.readline().split()
            if '$' in line or '!' in line:
                pass
            else:
                options.append(line)
        n_atoms = int(options[0][0])
        if int(options[0][1]) != 0:
            n_atoms += 1
        while True:
            line = lines.readline().split()
            if not line:
                break
            name = line[0]
            for i in range(n_atoms):
                line = lines.readline().split()
                if '!' in line:
                    pass
                else:
                    if len(line) == 4:
                        input_molecule[0].append(line[0])
                        input_molecule[1].append(line[1:])
                    elif len(line) == 5:
                        input_molecule[0].append(line[0])
                        input_molecule[1].append(line[1:-1])
                    else:
                        sys.exit('Wrong input format')
            if input_molecule:
                molecules.append(Molecule(input_molecule, name=name))
                input_molecule = [[], []]
    return [molecules, options]

# OUTPUT part
def write_wyfsym_measure(results, output_name):
    output = open('results/' + output_name + '.wout', 'w')

    results.print_alpha_mo_IRD()
    results.print_beta_mo_IRD()
    results.print_wf_mo_IRD()
    results.print_CSM()
    results.print_ideal_group_table()
    results.print_overlap_mo_alpha()
    results.print_overlap_mo_beta()
    results.print_overlap_wf()
    for i in range(results.dgroup):
        results.print_symmetry_operation_matrix(i)
        results.print_symmetry_transformed_coordinates(i)


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list
