from cosymlib.molecule import Molecule, Geometry, ElectronicStructure
from cosymlib.file_io.tools import extract_geometries
from cosymlib.tools import atomic_number_to_element

import numpy as np
import os, sys, re


def _non_blank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line


# Read INPUT files
def read_generic_structure_file(input_name, read_multiple=False):
    # print('Reading file {}...'.format(os.path.basename(input_name)))
    if os.stat(input_name).st_size == 0:
        raise FileExistsError('File {} is empty'.format(os.path.basename(input_name)))
    file_name, file_extension = os.path.splitext(input_name)
    if 'molden' in file_name:
        file_extension = ' molden'

    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get('get_molecule_from_file_' + file_extension[1:])
    if method is None:
        method = possibles.get('get_geometry_from_file_' + file_extension[1:])
    if not method:
        raise NotImplementedError("File not recognized")

    return method(input_name, read_multiple=read_multiple)


def get_geometry_from_file_xyz(file_name, read_multiple=False):
    """
    Reads a XYZ file and returns the geometry of all structures in it
    :param file_name: file name
    :param read_multiple: read multiple files if available
    :return: list of Geometry objects
    """
    input_molecule = [[], []]
    geometries = []
    with open(file_name, mode='r') as lines:
        lines.readline()
        name = lines.readline()
        if name.strip():
            name = name.split()[0]
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    float(line.split()[1])
                    input_molecule[0].append(line.split()[0])
                    input_molecule[1].append(line.split()[1:])
                except (ValueError, IndexError):
                    if input_molecule[0]:
                        geometries.append(Geometry(symbols=input_molecule[0],
                                                  positions=input_molecule[1],
                                                  name=name))
                    input_molecule = [[], []]
                    name = line.split()[0]

        molecule = Geometry(symbols=input_molecule[0],
                                  positions=input_molecule[1],
                                  name=name)
        if read_multiple:
            geometries.append(molecule)
        else:
            return molecule

    return geometries


def get_geometry_from_file_cor(file_name, read_multiple=False):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :param read_multiple: read multiple files if available
    :return: list of Geometry objects
    """
    input_molecule = [[], []]
    geometries = []
    with open(file_name, mode='r') as lines:
        line = lines.readline()
        name = line.replace('**FRAG**', '')[:-1]
        name = name.replace(' ', '')
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    float(line.split()[1])
                    if len(line.split()) == 4:
                        input_molecule[0].append(line.split()[0])
                        input_molecule[1].append(line.split()[1:])
                    elif len(line.split()) == 5:
                        input_molecule[0].append(line.split()[0])
                        input_molecule[1].append(line.split()[1:-1])
                    else:
                        float(line.split()[2])
                except (ValueError, IndexError):
                    if input_molecule:
                        geometries.append(Geometry(symbols=input_molecule[0],
                                                   positions=input_molecule[1],
                                                   name=name))
                    input_molecule = [[], []]
                    name = line.replace('**FRAG**', '')[:-1]
                    name = name.replace(' ', '')
        geometries.append(Geometry(symbols=input_molecule[0],
                                   positions=input_molecule[1],
                                   name=name))
        if not read_multiple:
            return geometries[0]
    return geometries


# TODO: This function should handle open and close shell, mo energies and more. A lot to improve!
def get_molecule_from_file_fchk(file_name, read_multiple=False):
    key_list = ['Charge', 'Multiplicity', 'Number of alpha electrons', 'Number of beta electrons', 'Atomic numbers',
                'Current cartesian coordinates', 'Shell type', 'Number of primitives per shell', 'Shell to atom map',
                'Primitive exponents', 'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha Orbital Energies', 'Alpha MO coefficients', 'Beta MO coefficients']
    input_molecule = [[] for _ in range(len(key_list))]
    read = False
    with open(file_name, mode='r') as lines:
        name = lines.readline().split()[0]#.strip('\n')
        line = lines.readline().split()
        basis_set = line[-1]
        if 'R' in line[1]:
            del key_list[-1]

        n = 1
        options = True
        for line in lines:
            if read:
                try:
                    float(line.split()[0])
                    input_molecule[n].append(line.split())
                except ValueError:
                    read = False

            for idn, key in enumerate(key_list):
                if key in line:
                    if n == len(key_list) - 1:
                        break
                    if options and idn != 4:
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

        for n in range(4, len(input_molecule)):
            input_molecule[n] = reformat_input(input_molecule[n])
        bohr_to_angstrom = 0.529177249
        coordinates = np.array(input_molecule[5], dtype=float).reshape(-1, 3) * bohr_to_angstrom

        atomic_number = [int(num) for num in input_molecule[4]]
        symbols = []
        for number in atomic_number:
            symbols.append(atomic_number_to_element(number))

        basis = basis_format(basis_set_name=basis_set,
                             atomic_numbers=atomic_number,
                             atomic_symbols=symbols,
                             shell_type=input_molecule[6],
                             n_primitives=input_molecule[7],
                             atom_map=input_molecule[8],
                             p_exponents=input_molecule[9],
                             c_coefficients=input_molecule[10],
                             p_c_coefficients=input_molecule[11])

        geometry = Geometry(symbols=symbols,
                            positions=coordinates,
                            name=name)

        Ca = np.array(input_molecule[13], dtype=float).reshape(-1, int(np.sqrt(len(input_molecule[13]))))
        energies_alpha = np.array(input_molecule[12], dtype=float).tolist()
        if input_molecule[14]:
            Cb = np.array(input_molecule[14], dtype=float).reshape(-1, int(np.sqrt(len(input_molecule[14]))))
        else:
            Cb = []

        electronic_structure = ElectronicStructure(charge=input_molecule[0][0],
                                                   multiplicity=input_molecule[1][0],
                                                   basis=basis,
                                                   orbital_coefficients=[Ca, Cb],
                                                   mo_energies=energies_alpha,
                                                   alpha_electrons=input_molecule[2],
                                                   beta_electrons=input_molecule[3])

        if read_multiple:
            return [Molecule(geometry, electronic_structure)]
        else:
            return Molecule(geometry, electronic_structure)


def get_molecule_from_file_molden(file_name, read_multiple=False):
    key_list = ['Charge', 'Multiplicity', 'Atomic numbers', 'Current cartesian coordinates',
                'Shell type', 'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha MO coefficients', 'Beta MO coefficients', 'MO Energies']

    type_list = {'s': '0',
                 'p': '1',
                 'd': '2',
                 'f': '3',
                 'sp': '-1'}

    input_molecule = {key: [] for key in key_list}
    read_molden = False
    read_coordinates = False
    read_basis = False
    read_coefficients = False
    occupation = {'Alpha': [], 'Beta': []}
    with open(file_name, mode='r') as lines:
        lines.readline()
        lines.readline()
        name = lines.readline()
        for line in _non_blank_lines(lines):

            if '[' in line:
                pass
            elif read_molden:
                if '======= END OF MOLDEN-FORMAT MOLECULAR ORBITALS =======' in line:
                    break
                if read_coordinates:
                    input_molecule['Atomic numbers'].append(line.split()[2])
                    input_molecule['Current cartesian coordinates'].append(line.split()[3:])
                if read_basis:
                    try:
                        number = float(line.split()[0].replace('D', 'E'))
                        if number-int(number) != 0.:
                            input_molecule['Primitive exponents'].append(line.split()[0].replace('D', 'E'))
                            input_molecule['Contraction coefficients'].append(line.split()[1].replace('D', 'E'))
                            if len(line.split()) > 2:
                                input_molecule['P(S=P) Contraction coefficients'].append(line.split()[2].replace('D', 'E'))
                            else:
                                input_molecule['P(S=P) Contraction coefficients'].append(0.)
                        else:
                            atom = line.split()[0]
                    except ValueError:
                        input_molecule['Shell to atom map'].append(atom)
                        input_molecule['Shell type'].append(type_list[line.split()[0].lower()])
                        input_molecule['Number of primitives per shell'].append(line.split()[1])
                if read_coefficients:
                    if 'Sym' in line:
                        pass
                    else:
                        if 'Ene' in line:
                            input_molecule['MO Energies'].append(float(re.split('[= ]', line)[-1]))
                        elif 'Spin' in line:
                            spin = re.split('[= ]', line)[-1]
                            input_molecule[spin + ' MO coefficients'].append([])
                        elif 'Occup' in line:
                            occupation[spin].append(float(line.split('=')[-1]))
                        else:
                            input_molecule[spin+' MO coefficients'][-1].append(line.split()[1])

            if '[Atoms]' in line:
                read_molden = True
                read_coordinates = True
                if 'AU' in line:
                    bohr_to_angstrom = 0.529177249
                else:
                    bohr_to_angstrom = 1
            elif '[GTO]' in line:
                read_coordinates = False
                read_basis = True
            elif '[MO]' in line:
                read_basis = False
                read_coefficients = True
            elif '[5D]' in line:
                input_molecule['Shell type'] = [x.replace('2', '-2') for x in input_molecule['Shell type']]
            elif '[7F]' in line:
                input_molecule['Shell type'] = [x.replace('3', '-3') for x in input_molecule['Shell type']]

        coordinates = np.array(input_molecule['Current cartesian coordinates'], dtype=float).reshape(-1, 3) * bohr_to_angstrom

        total_n_electrons = 0
        symbols = []
        for atom_number in input_molecule['Atomic numbers']:
            total_n_electrons += int(atom_number)
            symbols.append(tools.atomic_number_to_element(atom_number))

        occupation['Alpha'] = np.array(occupation['Alpha'], dtype=float)
        occupation['Beta'] = np.array(occupation['Beta'], dtype=float)
        if len(occupation['Beta']) == 0:
            if 2.0 not in occupation['Alpha']:
                occupation['Alpha'] = 2*occupation['Alpha']
            n_electrons = sum(occupation['Alpha'])
            input_molecule['Multiplicity'].append(sum(map(lambda x: x % 2 == 1, occupation['Alpha'])) + 1)
        else:
            if 2.0 not in occupation['Beta']:
                occupation['Beta'] = 2*occupation['Beta']
            n_electrons = sum(occupation['Alpha'] + occupation['Beta'])
            input_molecule['Multiplicity'].append(sum(occupation['Alpha'] - occupation['Beta']) + 1)

        input_molecule['Charge'].append(int(total_n_electrons - n_electrons))

        basis = basis_format(basis_set_name='UNKNOWN',
                             atomic_numbers=input_molecule['Atomic numbers'],
                             atomic_symbols=symbols,
                             shell_type=input_molecule['Shell type'],
                             n_primitives=input_molecule['Number of primitives per shell'],
                             atom_map=input_molecule['Shell to atom map'],
                             p_exponents=input_molecule['Primitive exponents'],
                             c_coefficients=input_molecule['Contraction coefficients'],
                             p_c_coefficients=input_molecule['P(S=P) Contraction coefficients'])

        geometry = Geometry(symbols=symbols,
                            positions=coordinates,
                            name=name)

        Ca = np.array(input_molecule['Alpha MO coefficients'], dtype=float)
        if input_molecule['Beta MO coefficients']:
            Cb = np.array(input_molecule['Beta MO coefficients'], dtype=float)
        else:
            Cb = []

        ee = ElectronicStructure(charge=input_molecule['Charge'][0],
                                 multiplicity=input_molecule['Multiplicity'][0],
                                 basis=basis,
                                 orbital_coefficients=[Ca, Cb],
                                 mo_energies=input_molecule['MO Energies'])

        if read_multiple:
            return [Molecule(geometry, ee)]
        else:
            Molecule(geometry, ee)


def get_geometry_from_file_ref(file_name, read_multiple=False):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :param read_multiple: read all geometries inside the file
    :return: list of Geometry objects
    """
    input_molecule = []
    structures = []
    with open(file_name, mode='r') as lines:
        lines.readline()
        lines.readline()
        name = lines.readline().split()[0]
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    if len(line.split()) > 3:
                        float(line.split()[1])
                        input_molecule.append(line.split()[1:])
                    else:
                        float(line.split()[0])
                        input_molecule.append(line.split())
                except (ValueError, IndexError):
                    if input_molecule:
                        structures.append(Geometry(positions=input_molecule,
                                                  name=name))
                    input_molecule = []
                    name = line.split()[0]
        structures.append(Geometry(positions=input_molecule,
                                   name=name))
    if read_multiple is False:
        return structures[0]

    return structures


def get_geometry_from_file_pdb(file_name, read_multiple=False):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :param read_multiple: read multiple files if available
    :return: list of Geometry objects
    """

    file_txt = open(file_name, mode='r').read()

    index_list = []
    for i in re.finditer('TITLE', file_txt):
        index_list.append(i.start())
    index_list.append(len(file_txt))

    geometries = []
    for i in range(len(index_list)-1):
        mol_section = file_txt[index_list[i]:index_list[i+1]]
        coordinates = []
        symbols = []
        connect = []
        name = ''
        for line in mol_section.split('\n'):
            if line.find('TITLE') > -1:
                name = ' '.join(line.split()[1:])
            if line.find('HETATM') > -1 or line.find('ATOM') > -1:
                coordinates.append([float(num) for num in line.split()[4:7]])
                symbols.append(line.split()[2])
            if line.find('CONECT') > -1:
                connect.append([int(num) for num in line.split()[1:]])
            if line.find('TER') > -1 or line.find('END') > -1:
                break

        connectivity = []
        for atom in connect:
            for j in atom[1:]:
                connectivity.append((atom[0], j))

        geometries.append(Geometry(symbols=symbols,
                                   positions=coordinates,
                                   name=name,
                                   connectivity=connectivity))

    if read_multiple is False:
        return geometries[0]

    return geometries


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):
    type_list = {'0': ['s', 1],
                 '1': ['p', 3],
                 '2': ['d', 6],
                 '3': ['f', 10],
                 '-1': ['sp', 4],
                 '-2': ['d_', 5],
                 '-3': ['f_', 7]}

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    atom_map = np.array(atom_map, dtype=int)

    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    # shell_type_index = [0] + np.cumsum([type_list['{}'.format(s)][1]
    #                                     for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = atomic_symbols[iatom]

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = type_list['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell + 1]

            shells_data.append({
                'shell_type': st[0],
                'p_exponents': p_exponents[ini_prim: fin_prim],
                'con_coefficients': c_coefficients[ini_prim: fin_prim],
                'p_con_coefficients': p_c_coefficients[ini_prim: fin_prim],
            })

        atoms_data.append({'shells': shells_data,
                           'symbol': symbol,
                           'atomic_number': atomic_number})

    basis_set['atoms'] = atoms_data

    return basis_set


def get_connectivity_from_file(file_name):

    file_txt = open(file_name, mode='r').read()

    connect = []
    for line in file_txt.split('\n'):
            connect.append([int(num) for num in line.split()])

    connectivity = []
    for atom in connect:
        for j in atom[1:]:
            connectivity.append((atom[0], j))

    return connectivity


# Get OUPUT files
def get_file_xyz_txt(structure):

    # Get xyz file format from Molecule or Geometry (or list)
    geometries = extract_geometries(structure, as_list=True)

    txt = ''
    for geometry in geometries:
        txt += '{}\n'.format(geometry.get_n_atoms())
        txt += '{}\n'.format(geometry.name if geometry.name is not None else '')
        for idp, position in enumerate(geometry.get_positions()):
            txt += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idp],
                                                                  position[0], position[1], position[2])
    return txt


# Support functions
def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list
