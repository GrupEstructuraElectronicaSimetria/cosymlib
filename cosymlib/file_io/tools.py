import yaml, sys
from cosymlib.molecule.geometry import Geometry
from cosymlib import __version__


def print_header(output=sys.stdout):
    output.write('-' * 70 + '\n')
    output.write(' COSYMLIB v{}\n Electronic Structure & Symmetry Group\n'.format(__version__))
    output.write(' Institut de Quimica Teorica i Computacional (IQTC)\n')
    output.write(' Universitat de Barcelona\n')
    output.write('-' * 70 + '\n\n')


def print_footer(output=sys.stdout):
    output.write('\n' + '-' * 70 + '\n')
    output.write(' ' * 20 + 'End of calculation\n')
    output.write('-' * 70 + '\n\n')


def print_input_info(initial_geometries, output=sys.stdout):

    for ids, geometry in enumerate(initial_geometries):
        output.write('Structure {} : {}\n'.format(ids+1, geometry.name))

        for idn, array in enumerate(geometry.get_positions()):
            output.write('{:2s}'.format(geometry.get_symbols()[idn]))
            output.write(' {:11.7f} {:11.7f} {:11.7f}'.format(array[0], array[1], array[2]))

            if geometry.get_connectivity() is not None:
                output.write('  ' +
                             '  '.join(['{:3}'.format(i+1) for i in range(geometry.get_n_atoms())
                                        if [idn+1, i+1] in geometry.get_connectivity()])+
                             '\n')
        output.write('\n')


def add_extra_keywords(args, filename):
    with open(filename, 'r') as stream:
        input_parameters = yaml.load(stream, Loader=yaml.FullLoader)
    try:
        for key, value in input_parameters.items():
            if key.lower() in args:
                setattr(args, key.lower(), value)
            else:
                raise KeyError("'{}' is not a valid input keyword".format(key))
    except AttributeError:
        raise KeyError('Incorrect input file format')


def extract_geometries(structure_list, as_list=False):
    """
    extracts the Geometry from Molecule or list of Molecule
    if Geometry or list of Geometry is introduced no changes are made

    :param structure_list:
    :param as_list:
    :return:
    """

    if isinstance(structure_list, list):
        geometry_list = []
        n_atoms_ref = None
        for structure in structure_list:
            if isinstance(structure, Geometry):
                geometry_list.append(structure)
                n_atoms = structure.get_n_atoms()
            else:
                geometry_list.append(structure.geometry)
                n_atoms = structure.geometry.get_n_atoms()

            if n_atoms_ref is not None:
                if n_atoms != n_atoms_ref:
                    raise Exception('All geometries should contain the same number of atoms!')
                n_atoms_ref = structure.get_n_atoms()

        return geometry_list

    else:
        if isinstance(structure_list, Geometry):
            return [structure_list] if as_list else structure_list
        else:
            return [structure_list.geometry] if as_list else structure_list.geometry

