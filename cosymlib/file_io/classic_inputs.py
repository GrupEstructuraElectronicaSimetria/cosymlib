import sys
from cosymlib.molecule import Geometry
from cosymlib.file_io import read_generic_structure_file
import os.path


def read_old_input(file_name):
    """
    Reads the old Shape's program input
    :param file_name: file name
    :return: list of Geometry objects and options
    """
    input_molecule = [[], []]
    structures = []
    options = {'%out': None, '%conquest': None, '%external': False, '%fullout': False, '%test': False,
               '%n_atoms': 0, '%central_atom': 0, '%labels': 0, '%path': False}

    with open(file_name, mode='r') as lines:
        while True:
            line = lines.readline().split()
            if '$' in line or '!' in line:
                pass
            elif any('%' in word for word in line):
                if len(line) > 1:
                    options[line[0]] = line[1]
                else:

                    options[line[0]] = True
            else:
                try:
                    int(line[0])
                    if options['%n_atoms'] == 0:
                        options['%n_atoms'] = int(line[0])
                        options['%central_atom'] = int(line[1])
                    else:
                        options['%labels'] = line
                except (ValueError, IndexError):
                    break

        n_atoms = options['%n_atoms']
        if options['%central_atom'] != 0:
            n_atoms += 1
        if options['%conquest'] is not None:
            dir = os.path.dirname(os.path.abspath(file_name))
            structures = read_generic_structure_file(os.path.join(dir, options['%conquest'] + '.cor'), read_multiple=True)
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
