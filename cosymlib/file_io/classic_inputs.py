from cosymlib.file_io import get_geometry_from_file_cor
from cosymlib.file_io import errors

import os
import tempfile
import warnings


def read_old_input(file_name):
    """
    Reads the old Shape's program input
    :param file_name: file name
    :return: list of Geometry objects and options
    """

    options = {'%out': None, '%conquest': None, '%external': False, '%fullout': False, '%test': False,
               '%n_atoms': 0, '%central_atom': 0, '%labels': 0, '%path': False}

    idl = 0
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
            idl += 1

        n_atoms = options['%n_atoms']
        if options['%central_atom'] != 0:
            n_atoms += 1
        if options['%conquest'] is not None:
            dir = os.path.dirname(os.path.abspath(file_name))
            structures = get_geometry_from_file_cor(os.path.join(dir, options['%conquest'] + '.cor'), read_multiple=True)
        else:
            tmp = tempfile.NamedTemporaryFile(mode='w+t', dir=os.getcwd())
            tmp_lines = lines.readlines()
            try:
                # Write data to the temporary file
                tmp.write(line[0]+'\n')
                idl += 1
                for line in tmp_lines:
                    if line.strip() == '':
                        warnings.warn('Line {} is empty'.format(idl + 1), errors.EmptyLineWarning)
                    else:
                        tmp.write(line)
                    idl += 1
                tmp.seek(0)
                structures = get_geometry_from_file_cor(tmp.name, read_multiple=True)
            finally:
                tmp.close()

    return [structures, options]
