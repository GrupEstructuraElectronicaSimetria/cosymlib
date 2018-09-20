import os
import sys
from symeess.molecule import Molecule, Geometry
import numpy as np


# INPUT part
def read_input_file(input_name):

    file_name, file_extension = os.path.splitext(input_name)
    method_name = 'read_file_' + file_extension[1:]
    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(method_name)
    if not method:
        raise NotImplementedError("Method %s not implemented" % method_name)

    return method(input_name)


def read_geometry_from_xyz_file(file_name):
    """
    Reads a XYZ file and returns the geometry of all structures in it
    :param file_name: file name
    :return: list of Geometry objects
    """
    input_molecule = [[], []]
    structures = []
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
                    if input_molecule:
                        structures.append(Geometry(input_molecule, name=name))
                    input_molecule = [[], []]
                    name = lines.readline().split()[0]
        structures.append(Geometry(input_molecule, name=name))
    return structures


def read_geometry_from_cor_file(file_name):
    """
    Reads a Conquest formatted file and the geometry of all structures in it
    :param file_name: file name
    :return: list of Geometry objects
    """
    input_molecule = [[], []]
    structures = []
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
                        structures.append(Geometry(input_molecule, name=name))
                    input_molecule = [[], []]
                    name = line.split()[0]
        structures.append(Geometry(input_molecule, name=name))
    return structures


def read_old_input(file_name):
    """
    Reads the old Shape's program input
    :param file_name: file name
    :return: list of Geometry objects and
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
                molecules.append(Molecule(input_molecule, name=name))
                input_molecule = [[], []]
    return [molecules, options]


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
        return Molecule(structure=input_molecule[3:5],
                        ee=input_molecule[:3] + input_molecule[5:],
                        name=name)


# OUTPUT part
def write_wfnsym_measure(label, geometry, wfnsym_results, output_name):

    if not os.path.exists('./results'):
        os.makedirs('./results')
    output = open('results/' + output_name + '.wout', 'w')

    # Print Outputs
    output.write('MEASURES OF THE SYMMETRY GROUP:   {}\n'.format(label))
    output.write('Basis: {}\n'.format('6-31G(d)'))
    output.write('--------------------------------------------\n')
    output.write(' Atomic Coordinates (Angstroms)\n')
    output.write('--------------------------------------------\n')
    for array in geometry.get_positions():
        output.write(' {:10.7f} {:10.7f} {:10.7f}\n'.format(array[0], array[1], array[2]))
    output.write('--------------------------------------------\n')
    for i in range(wfnsym_results.dgroup):
        output.write('\n')
        output.write('@@@ Operation {0}: {1}'.format(i + 1, wfnsym_results.SymLab[i]))
        output.write('\nSymmetry Transformation matrix\n')
        for array in wfnsym_results.SymMat[i]:
            output.write(' {:10.7f} {:10.7f} {:10.7f}\n'.format(array[0], array[1], array[2]))
        output.write('\n')
        output.write('Symmetry Transformed Atomic Coordinates (Angstroms)\n')
        for array in np.dot(geometry.get_positions(), wfnsym_results._SymMat[i].T):
            output.write(' {:10.7f} {:10.7f} {:10.7f}\n'.format(array[0], array[1], array[2]))

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


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list
