import os
import sys
from symeess.molecule import Molecule, Geometry, ElectronicStructure
import numpy as np
from symeess import tools


# INPUT part
def read_input_file(input_name):
    print('Reading file {}...'.format(os.path.basename(input_name)))
    if os.stat(input_name).st_size == 0:
        raise FileExistsError('File {} is empty'.format(os.path.basename(input_name)))
    file_name, file_extension = os.path.splitext(input_name)
    method_name = 'get_molecule_from_file_' + file_extension[1:]
    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(method_name)
    if not method:
        raise NotImplementedError("Method %s not implemented" % method_name)

    return method(input_name)


def get_molecule_from_file_xyz(file_name):
    """
    Reads a XYZ file and returns the geometry of all structures in it
    :param file_name: file name
    :return: list of Geometry objects
    """
    input_molecule = [[], []]
    molecules = []
    with open(file_name, mode='r') as lines:
        lines.readline()
        name = lines.readline().split()[0]
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
                        molecules.append(Geometry(symbols=input_molecule[0],
                                                  positions=input_molecule[1],
                                                  name=name))
                    input_molecule = [[], []]
                    name = line.split()[0]

        molecules.append(Geometry(symbols=input_molecule[0],
                                  positions=input_molecule[1],
                                  name=name))

    return molecules


def get_molecule_from_file_cor(file_name):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :return: list of Geometry objects
    """
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
                    if len(line.split()) == 4:
                        input_molecule[0].append(line.split()[0])
                        input_molecule[1].append(line.split()[1:])
                    elif len(line.split()) == 5:
                        input_molecule[0].append(line.split()[0])
                        input_molecule[1].append(line.split()[1:-1])
                    else:
                        sys.exit('Wrong input format')
                    # input_molecule[0].append(line.split()[0])
                    # input_molecule[1].append(line.split()[1:-1])
                except (ValueError, IndexError):
                    if input_molecule:
                        molecules.append(Geometry(symbols=input_molecule[0],
                                                  positions=input_molecule[1],
                                                  name=name))
                    input_molecule = [[], []]
                    name = line.split()[0]
        molecules.append(Geometry(symbols=input_molecule[0],
                                  positions=input_molecule[1],
                                  name=name))
    return molecules


def get_molecule_from_file_fchk(file_name):
    key_list = ['Charge', 'Multiplicity', 'Atomic numbers', 'Current cartesian coordinates',
                'Shell type', 'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha MO coefficients', 'Beta MO coefficients']
    input_molecule = [[] for _ in range(len(key_list))]
    read = False
    with open(file_name, mode='r') as lines:
        name = lines.readline()
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
                    # input_molecule[n] = reformat_input(input_molecule[n])
                    read = False

            for idn, key in enumerate(key_list):
                if key in line:
                    if n == len(key_list) - 1:
                        break
                    if options and idn != 2:
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

        for n in range(2, len(input_molecule) - 1):
            input_molecule[n] = reformat_input(input_molecule[n])
        bohr_to_angstrom = 0.529177249
        coordinates = np.array(input_molecule[3], dtype=float).reshape(-1, 3) * bohr_to_angstrom

        atomic_number = [int(num) for num in input_molecule[2]]
        symbols = []
        for number in atomic_number:
            symbols.append(tools.atomic_number_to_element(number))

        basis = basis_format(basis_set_name=basis_set,
                             atomic_numbers=atomic_number,
                             atomic_symbols=symbols,
                             shell_type=input_molecule[4],
                             n_primitives=input_molecule[5],
                             atom_map=input_molecule[6],
                             p_exponents=input_molecule[7],
                             c_coefficients=input_molecule[8],
                             p_c_coefficients=input_molecule[9])

        geometry = Geometry(symbols=input_molecule[2],
                            positions=coordinates,
                            name=name)

        ee = ElectronicStructure(charge=input_molecule[0][0],
                                 multiplicity=input_molecule[1][0],
                                 basis=basis,
                                 orbital_coefficients=[input_molecule[10], input_molecule[11]])

        return [Molecule(geometry, ee)]

def get_molecule_from_file_ref(file_name):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :return: list of Geometry objects
    """
    input_molecule = []
    molecules = []
    with open(file_name, mode='r') as lines:
        lines.readline()
        lines.readline()
        name = lines.readline().split()[0]
        for line in lines:
            if '$' in line or '#' in line:
                pass
            else:
                try:
                    float(line.split()[0])
                    input_molecule.append(line.split())
                except (ValueError, IndexError):
                    if input_molecule:
                        molecules.append(Geometry(positions=input_molecule,
                                                  name=name))
                    input_molecule = []
                    name = line.split()[0]
        molecules.append(Geometry(positions=input_molecule,
                                  name=name))
    return molecules


def read_old_input(file_name):
    """
    Reads the old Shape's program input
    :param file_name: file name
    :return: list of Geometry objects and options
    """
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
            if input_molecule[0]:
                molecules.append(Geometry(symbols=input_molecule[0],
                                          positions=input_molecule[1],
                                          name=name))
                input_molecule = [[], []]
    return [molecules, options]


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):
    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d', 5],
                '-3': ['f', 7]}

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    atom_map = np.array(atom_map, dtype=int)

    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    shell_type_index = [0] + np.cumsum([typeList['{}'.format(s)][1]
                                        for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = atomic_symbols[iatom]

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
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


# OUTPUT part
def header(output):
    output.write('-' * 70 + '\n')
    output.write('SYMEESS v0.6.3 \n'
                 'Electronic Structure Group,  Universitat de Barcelona\n')
    output.write('-' * 70 + '\n' + '\n')


def write_input_info(initial_geometries, output_name=None):
    if output_name is not None:
        # if not os.path.exists('./results'):
        #     os.makedirs('./results')
        output = open(output_name + '.tab', 'w')
    else:
        output = sys.stdout
    header(output)

    for ids, geometry in enumerate(initial_geometries):
        output.write('Structure {} : {}\n'.format(ids+1, geometry.get_name()))

        for idn, array in enumerate(geometry.get_positions()):
            output.write('{:2s}'.format(geometry.get_symbols()[idn]))
            output.write(' {:11.7f} {:11.7f} {:11.7f}\n'.format(array[0], array[1], array[2]))
        output.write('\n')


def write_symgroup_measure(label, geometries, symgroup_results, output_name):
    if output_name is not None:
        # if not os.path.exists('./results'):
        #     os.makedirs('./results')
        output = open(output_name + '.zout', 'w')
        output2 = open(output_name + '.ztab', 'w')
    else:
        output = sys.stdout
        output2 = sys.stdout
    header(output)

    output.write('Evaluating symmetry operation : {}\n'.format(label))

    for idx, geometry in enumerate(geometries):
        output.write('{}\n'.format(geometry.get_name()))
        output.write('\n')
        output.write('Centered Structure\n')
        output.write('--------------------------------------------\n')
        center_mass = tools.center_mass(geometry.get_symbols(), geometry.get_positions())
        for idn, array in enumerate(geometry.get_positions()):
            array = array - center_mass
            output.write('{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(geometry.get_symbols()[idn],
                                                                    array[0], array[1], array[2]))
        output.write('--------------------------------------------\n')

        output.write('Optimal permutation\n')
        for idn, permutation in enumerate(symgroup_results[idx].optimum_permutation):
            output.write('{:2} {:2}\n'.format(idn + 1, permutation))
        output.write('\n')

        output.write('Inverted structure\n')
        for idn, axis in enumerate(symgroup_results[idx].nearest_structure):
            output.write('{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(geometry.get_symbols()[idn],
                                                                    axis[0], axis[1], axis[2]))
        output.write('\n')

        output.write('Reference axis\n')
        for array in symgroup_results[idx].reference_axis:
            output.write('{:12.8f} {:12.8f} {:12.8f}\n'.format(array[0], array[1], array[2]))
        output.write('\n')

        output.write('Symmetry measure {:.5f}\n'.format(symgroup_results[idx].csm))
        output.write('..................................................\n')
    output2.write('Evaluating symmetry operation : {}\n \n'.format(label))
    for idx, geometry in enumerate(geometries):
        output2.write('{:>5} {:10.5f}\n'.format(geometry.get_name(), symgroup_results[idx].csm))

    output.close()
    output2.close()


def write_wfnsym_measure(label, molecule, wfnsym_results, output_name):
    if output_name is not None:
        # if not os.path.exists('./results'):
        #     os.makedirs('./results')
        output = open(output_name + '.wout', 'w')
    else:
        output = sys.stdout
    header(output)

    geometry = molecule.geometry
    output.write('MEASURES OF THE SYMMETRY GROUP:   {}\n'.format(label))
    output.write('Basis: {}\n'.format(list(molecule.electronic_structure.basis.keys())[0]))
    output.write('--------------------------------------------\n')
    output.write(' Atomic Coordinates (Angstroms)\n')
    output.write('--------------------------------------------\n')
    for idn, array in enumerate(geometry.get_positions()):
        output.write('{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                                array[0], array[1], array[2]))
    output.write('--------------------------------------------\n')
    for i, label in enumerate(wfnsym_results.SymLab):
        output.write('\n')
        output.write('@@@ Operation {0}: {1}'.format(i + 1, wfnsym_results.SymLab[i]))
        output.write('\nSymmetry Transformation matrix\n')
        for array in wfnsym_results.SymMat[i]:
            output.write(' {:11.6f} {:11.6f} {:11.6f}\n'.format(array[0], array[1], array[2]))
        output.write('\n')
        output.write('Symmetry Transformed Atomic Coordinates (Angstroms)\n')

        for idn, array in enumerate(geometry.get_positions()):
            array2 = np.dot(array, wfnsym_results._SymMat[i].T)
            output.write('{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                                    array2[0], array2[1], array2[2]))

    output.write('\nIdeal Group Table\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.SymLab]))
    output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    for i, line in enumerate(wfnsym_results.ideal_gt):
        output.write('{:4}'.format(wfnsym_results.IRLab[i]) + '  '.join(['{:7.3f}'.format(s) for s in line]))
        output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')

    output.write('\nAlpha MOs: Symmetry Overlap Expectation Values\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.SymLab]))
    output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')

    for i, line in enumerate(wfnsym_results.mo_SOEVs_a):
        output.write('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))
        output.write('\n')

    output.write('\nBeta MOs: Symmetry Overlap Expectation Values\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.SymLab]))
    output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    for i, line in enumerate(wfnsym_results.mo_SOEVs_b):
        output.write('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))
        output.write('\n')

    output.write('\nWaveFunction: Symmetry Overlap Expectation Values\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.SymLab]))
    output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('a-wf' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_SOEVs_a]))
    output.write('\n')
    output.write('b-wf' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_SOEVs_b]))
    output.write('\n')
    output.write('WFN ' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_SOEVs]))
    output.write('\n')

    output.write('\nWaveFunction: CSM-like values\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.SymLab]))
    output.write('\n')
    output.write('   -------------------------------------------------------------------------------'
                 '------------------------------------------------------------------------\n')

    output.write('Grim' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.grim_coef]))
    output.write('\n')
    output.write('CSM ' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.csm_coef]))
    output.write('\n')

    output.write('\nAlpha MOs: Irred. Rep. Decomposition\n')
    output.write('   ---------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.IRLab]))
    output.write('\n')
    output.write('   ---------------------------------------------\n')
    for i, line in enumerate(wfnsym_results.mo_IRd_a):
        output.write('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))
        output.write('\n')

    output.write('\nBeta MOs: Irred. Rep. Decomposition\n')
    output.write('   ---------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.IRLab]))
    output.write('\n')
    output.write('   ---------------------------------------------\n')
    for i, line in enumerate(wfnsym_results.mo_IRd_b):
        output.write('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))
        output.write('\n')

    output.write('\nWaveFunction: Irred. Rep. Decomposition\n')
    output.write('   ---------------------------------------------\n')
    output.write('     ' + '  '.join(['{:^7}'.format(s) for s in wfnsym_results.IRLab]))
    output.write('\n')
    output.write('   ---------------------------------------------\n')
    output.write('a-wf' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_IRd_a]))
    output.write('\n')
    output.write('b-wf' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_IRd_b]))
    output.write('\n')
    output.write('WFN ' + '  '.join(['{:7.3f}'.format(s) for s in wfnsym_results.wf_IRd]))
    output.write('\n')
    print('WFNSYM : Calculation has finished normally ')
    output.close()


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list


# def bhor2a(array):
#     new_array = []
#     for xyz in array:
#         new_array.append(float(xyz) / 1.889726124993590)
#     return np.array(new_array)
