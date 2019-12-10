import sys
from cosym.molecule import Geometry
from cosym.file_io import read_input_file
import os.path


def read_old_input(file_name):
    """
    Reads the old Shape's program input
    :param file_name: file name
    :return: list of Geometry objects and options
    """
    input_molecule = [[], []]
    structures = []
    options = {'%out': None, '%conquest': None, '%external': False, '%fullout': False, '%reference': False,
               '%test': False, '%n_atoms': 0, '%central_atom': 0, '%labels': 0}

    with open(file_name, mode='r') as lines:
        while True:
            line = lines.readline()
            if '$' in line or '!' in line:
                pass
            elif '%' in line:
                if len(line) > 0:
                    options[line.split()[0]] = line.split()[1]
                else:
                    options[line.split()[0]] = True
            else:
                try:
                    int(line.split()[0])
                    if options['%n_atoms'] == 0:
                        options['%n_atoms'] = int(line.split()[0])
                        options['%central_atom'] = int(line.split()[1])
                    else:
                        options['%labels'] = line.split()
                except (ValueError, IndexError):
                    break

        n_atoms = options['%n_atoms']
        if options['%central_atom'] != 0:
            n_atoms += 1
        if options['%conquest'] is not None:
            dir = os.path.dirname(file_name)+'/'
            structures = read_input_file(dir+options['%conquest']+'.cor')
        else:
            while True:
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
                if input_molecule[0]:
                    structures.append(Geometry(symbols=input_molecule[0],
                                               positions=input_molecule[1],
                                               name=name))
                    input_molecule = [[], []]
                line = lines.readline().split()

    return [structures, options]
