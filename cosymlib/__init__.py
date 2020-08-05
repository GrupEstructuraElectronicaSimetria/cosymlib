__version__ = '0.8.5'

from cosymlib.molecule import Molecule, Geometry
from cosymlib import file_io
from cosymlib import tools
from cosymlib.utils import get_shape_path, plot_molecular_orbital_diagram, plot_symmetry_energy_evolution
from cosymlib.shape.tools import get_structure_references
import matplotlib.pyplot as plt

import sys
import numpy as np


def _get_symgroup_arguments(locals):
    kwargs = dict(locals)
    del kwargs['self']
    for element in ['self', 'output']:
        if element in kwargs:
            del kwargs[element]

    return kwargs


class Cosymlib:
    """
    Main class of cosymlib program that can perform all the jobs
    """

    def __init__(self,
                 structures,
                 ignore_atoms_labels=False,
                 ignore_connectivity=False,
                 connectivity=None,
                 connectivity_thresh=None):

        self._molecules = []
        if isinstance(structures, list):
            for structure in structures:
                if isinstance(structure, Molecule):
                    self._molecules.append(structure)
                elif isinstance(structure, Geometry):
                    self._molecules.append(Molecule(structure))
                else:
                    raise AttributeError('Molecule object not found')
        else:
            if isinstance(structures, Molecule):
                self._molecules.append(structures)
            elif isinstance(structures, Geometry):
                self._molecules.append(Molecule(structures))
            else:
                raise AttributeError('Molecule object not found')

        for molecule in self._molecules:
            if ignore_atoms_labels:
                molecule.geometry.set_symbols('X' * molecule.geometry.get_n_atoms())
            if connectivity_thresh is not None:
                molecule.geometry.generate_connectivity(thresh=connectivity_thresh)
            if connectivity is not None:
                molecule.geometry.set_connectivity(connectivity)
            if ignore_connectivity:
                molecule.geometry.set_connectivity(None)

    def get_n_atoms(self):
        n_atoms_unique_list = np.unique([mol.geometry.get_n_atoms() for mol in self._molecules])
        if len(n_atoms_unique_list) > 1:
            raise Exception('Not all structures have same number of atoms')

        return n_atoms_unique_list[0]

    def get_geometries(self):
        return [mol.geometry for mol in self._molecules]

    def print_info(self):
        print('\033[1m{:20}   {:^5}\033[0m'.format('name', 'atoms'))
        for molecule in self._molecules:
            print('{:20} : {:5}'.format(molecule.name, molecule.geometry.get_n_atoms()))
        print('Total structures: {}'.format(len(self._molecules)))

    @property
    def molecules(self):
        return self._molecules

    def print_shape_measure(self, shape_reference, central_atom=0, fix_permutation=False, output=sys.stdout):
        """
        Method that prints to file shape's measure


        :param shape_reference: reference label, list of labels, 'all' or polyhedra list
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's measure in the output_name.tab file
        """

        if shape_reference == 'all':
            vertices = self.get_n_atoms() - int(bool(central_atom))
            reference_list = get_structure_references(vertices)
        else:
            if isinstance(shape_reference, (str, Geometry)):
                reference_list = [shape_reference]
            else:
                reference_list = shape_reference

        molecules_names = [molecule.name for molecule in self._molecules]

        measure_list = []
        references_names = []
        for reference in reference_list:
            measure_list.append(self.get_shape_measure(reference,
                                                       'measure',
                                                       central_atom,
                                                       fix_permutation))

            if type(reference) is Geometry:
                references_names.append(reference.name)
            else:
                references_names.append(reference)

        txt_shape = '{}'.format('Structure')
        max_name = len(max(molecules_names, key=len))
        if max_name < 9:
            n = 5
        else:
            n = max_name - 4
        for label in references_names:
            n += len(label)
            txt_shape += '{}'.format(label.rjust(n))
            n = 12 - len(label)
        txt_shape += '\n\n'

        for idx, name in enumerate(molecules_names):
            max_name = len(max(molecules_names, key=len))
            txt_shape += '{}'.format(name)
            if max_name < 9:
                n = 18 - len(name)
            else:
                n = 9 + max_name - len(name)
            for idn, label in enumerate(references_names):
                txt_shape += ', {:{width}.{prec}f}'.format(measure_list[idn][idx], width=n, prec=3)
                n = 11
            txt_shape += '\n'
        txt_shape += '\n'

        output.write(txt_shape)

    def print_shape_structure(self, shape_reference, central_atom=0, fix_permutation=False, output=sys.stdout):
        """
        Method that prints to file shape's structure

        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                                Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """

        if shape_reference == 'all':
            vertices = self.get_n_atoms() - int(bool(central_atom))
            reference_list = get_structure_references(vertices)
        else:
            if isinstance(shape_reference, (str, Geometry)):
                reference_list = [shape_reference]
            else:
                reference_list = shape_reference

        shape_results_structures = []
        references = []
        for reference in reference_list:
            if isinstance(reference, str):
                references.append(reference)

            else:
                references.append(reference.name)
                reference = reference.get_positions()

            shape_results_structures.append(self.get_shape_measure(reference,
                                                                   'structure',
                                                                   central_atom,
                                                                   fix_permutation))
        geometries = []
        for idm, molecule in enumerate(self._molecules):
            geometries.append(molecule.geometry)
            for idl, reference in enumerate(references):
                 shape_results_structures[idl][idm].set_name(molecule.name + '_' + reference)
                 geometries.append(shape_results_structures[idl][idm])

        for geometry in geometries:
            output.write(file_io.get_file_xyz_txt(geometry))

    # TODO: Change name of all symgroup named functions
    def print_geometric_measure_info(self, label, multi=1, central_atom=0, center=None, output=sys.stdout):
        kwargs = _get_symgroup_arguments(locals())
        sep_line = '..................................................\n'

        txt = 'Evaluating symmetry operation : {}\n'.format(label)

        for idx, molecule in enumerate(self._molecules):
            molecule.geometry._symmetry.set_parameters(kwargs)
            txt += '{}\n'.format(molecule.name)
            txt += '\n'
            txt += 'Centered Structure\n'
            txt += sep_line
            center_mass = tools.center_mass(molecule.geometry.get_symbols(), molecule.geometry.get_positions())
            for idn, array in enumerate(molecule.geometry.get_positions()):
                array = array - center_mass
                txt += '{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(molecule.geometry.get_symbols()[idn],
                                                                      array[0], array[1], array[2])
            txt += sep_line

            txt += 'Optimal permutation\n'
            for idn, permutation in enumerate(molecule.geometry._symmetry.optimum_permutation(label)):
                txt += '{:2} {:2}\n'.format(idn + 1, permutation)
            txt += '\n'

            txt += 'Inverted structure\n'
            for idn, axis in enumerate(molecule.geometry._symmetry.nearest_structure(label)):
                txt += '{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(molecule.geometry.get_symbols()[idn],
                                                                      axis[0], axis[1], axis[2])
            txt += '\n'

            txt += 'Reference axis\n'
            for array in molecule.geometry._symmetry.reference_axis(label):
                txt += '{:12.8f} {:12.8f} {:12.8f}\n'.format(array[0], array[1], array[2])
            txt += '\n'
            txt += 'Symmetry measure {:.5f}\n'.format(molecule.geometry.get_symmetry_measure(**kwargs))
            txt += sep_line

        output.write(txt)

    def print_geometric_symmetry_measure(self, label, central_atom=0, center=None,
                                         output=sys.stdout):
        kwargs = _get_symgroup_arguments(locals())

        txt = 'Evaluating symmetry operation : {}\n \n'.format(label)
        for idx, molecule in enumerate(self._molecules):
            csm = molecule.geometry.get_symmetry_measure(**kwargs)
            max_name = len(max(molecule.name, key=len))
            txt += '{} '.format(molecule.name)
            if max_name < 9:
                n = 18 - len(molecule.name)
            else:
                n = 9 + max_name - len(molecule.name)
            txt += '{:{width}.{prec}f}\n'.format(csm, width=n, prec=3)

        output.write(txt)

    def print_symmetry_nearest_structure(self, label, central_atom=0, center=None,
                                         output=sys.stdout):
        kwargs = _get_symgroup_arguments(locals())

        for idm, molecule in enumerate(self._molecules):
            geometry = molecule.geometry.get_symmetry_nearest_structure(**kwargs)
            output.write(file_io.get_file_xyz_txt(geometry))

    # This should be substituted by calling methods within this class
    def OLD_print_wnfsym_measure_verbose(self, group, axis=None, axis2=None, center=None, output=sys.stdout):
        self.print_wnfsym_sym_matrices(group, axis=axis, axis2=axis2, center=center, output=output)
        self.print_wnfsym_irreducible_repr(group, axis=axis, axis2=axis2, center=center, output=output)

    def print_electronic_symmetry_measure(self, group, axis=None, axis2=None, center=None, output=sys.stdout):

        txt = ''
        first = True
        for molecule in self._molecules:
            wf_measure = molecule.get_wf_symmetry(group, axis=axis, axis2=axis2, center=center)

            if first:
                sep_line = '          ' + '---------' * len(wf_measure['labels']) + '\n'

                txt += '\nWaveFunction: CSM-like values\n'
                txt += sep_line
                txt += '           ' + '  '.join(['{:^7}'.format(s) for s in wf_measure['labels']])
                txt += '\n'
                txt += sep_line

            txt += '{:<9} '.format(molecule.name) + '  '.join(['{:7.3f}'.format(s) for s in wf_measure['csm']])
            txt += '\n'
            first = False

        output.write(txt)

    def print_electronic_density_measure(self, group, axis=None, axis2=None, center=None, output=sys.stdout):

        txt = ''
        txt2 = ''
        first = True
        for molecule in self._molecules:
            dens_measure = molecule.get_dens_symmetry(group, axis=axis, axis2=axis2, center=center)

            if first:
                sep_line = '          ' + '---------' * len(dens_measure['labels']) + '\n'

                txt += '\nDensity: CSM-like values\n'
                txt += sep_line
                txt += '           ' + '  '.join(['{:^7}'.format(s) for s in dens_measure['labels']])
                txt += '\n'
                txt += sep_line

                txt2 += '--------------\n'
                txt2 += 'Total CSM {}\n'.format(group)
                txt2 += '--------------\n'

            txt += '{:<9} '.format(molecule.name) + '  '.join(['{:7.3f}'.format(s) for s in dens_measure['csm_coef']])
            txt += '\n'
            first = False
            axes_information = molecule.get_symmetry_axes(group, axis=axis, axis2=axis2, center=center)
            txt2 += '{:<9} '.format(molecule.name) + '{:7.3f}\n'.format(dens_measure['csm'])
            txt2 += '\ncenter: ' + '  '.join(['{:12.8f}'.format(s) for s in axes_information['center']])
            txt2 += '\n'
            txt2 += 'axis  : ' + '  '.join(['{:12.8f}'.format(s) for s in axes_information['axis']])
            txt2 += '\n'
            txt2 += 'axis2 : ' + '  '.join(['{:12.8f}'.format(s) for s in axes_information['axis2']])

        txt += '\n' + txt2

        output.write(txt)

    def print_wnfsym_sym_matrices(self, group, axis=None, axis2=None, center=None, output=sys.stdout):

        txt = ''
        for molecule in self._molecules:

            sym_mat = molecule.get_symmetry_matrix(group, axis=axis, axis2=axis2, center=center)

            geometry = molecule.geometry
            sym_lables = sym_mat['labels']
            sym_matrices = sym_mat['matrix']

            sep_line = '--------------------------------------------\n'
            txt += 'MEASURES OF THE SYMMETRY GROUP:   {}\n'.format(group)
            txt += 'Basis: {}\n'.format(list(molecule.electronic_structure.basis.keys())[0])
            txt += sep_line
            txt += ' Atomic Coordinates (Angstroms)\n'
            txt += sep_line
            for idn, array in enumerate(geometry.get_positions()):
                txt += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                                      array[0], array[1], array[2])
            txt += sep_line
            for i, group in enumerate(sym_lables):
                txt += '\n'
                txt += '@@@ Operation {0}: {1}'.format(i + 1, group)
                txt += '\nSymmetry Transformation matrix\n'
                for array in sym_matrices[i]:
                    txt += ' {:11.6f} {:11.6f} {:11.6f}\n'.format(array[0], array[1], array[2])
                txt += '\n'
                txt += 'Symmetry Transformed Atomic Coordinates (Angstroms)\n'

                for idn, array in enumerate(geometry.get_positions()):
                    array2 = np.dot(array, sym_matrices[i].T)
                    txt += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                                          array2[0], array2[1], array2[2])
        output.write(txt)

    def print_wnfsym_sym_ovelap(self, group, axis=None, axis2=None, center=None, output=sys.stdout):

        txt = ''
        for molecule in self._molecules:
            mo_overlap = molecule.get_mo_overlaps(group, axis=axis, axis2=axis2, center=center)
            wf_overlap = molecule.get_wf_overlaps(group, axis=axis, axis2=axis2, center=center)

            ideal_gt = molecule.get_ideal_group_table(group, axis=axis, axis2=axis2, center=center)
            wf_measure = molecule.get_wf_symmetry(group, axis=axis, axis2=axis2, center=center)

            sep_line = '     ' + '---------' * len(ideal_gt['ir_labels']) + '\n'

            txt += '\nMolecule : {}\n'.format(molecule.name)

            txt += '\nIdeal Group Table\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ideal_gt['labels']])
            txt += '\n'
            txt += sep_line
            for i, line in enumerate(ideal_gt['table']):
                txt += '{:4}'.format(ideal_gt['ir_labels'][i]) + '  '.join(['{:7.3f}'.format(s) for s in line])
                txt += '\n'
            txt += sep_line

            txt += '\nAlpha MOs: Symmetry Overlap Expectation Values\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ideal_gt['labels']])
            txt += '\n'
            txt += sep_line

            for i, line in enumerate(mo_overlap['alpha']):
                txt += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
                txt += '\n'

            txt += '\nBeta MOs: Symmetry Overlap Expectation Values\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ideal_gt['labels']])
            txt += '\n'
            txt += sep_line
            for i, line in enumerate(mo_overlap['beta']):
                txt += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
                txt += '\n'

            txt += '\nWaveFunction: Symmetry Overlap Expectation Values\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ideal_gt['labels']])
            txt += '\n'
            txt += sep_line
            txt += 'a-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_overlap['alpha']])
            txt += '\n'
            txt += 'b-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_overlap['beta']])
            txt += '\n'
            txt += 'WFN ' + '  '.join(['{:7.3f}'.format(s) for s in wf_overlap['total']])
            txt += '\n'

            txt += '\nWaveFunction: CSM-like values\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in wf_measure['labels']])
            txt += '\n'
            txt += sep_line

            txt += 'Grim' + '  '.join(['{:7.3f}'.format(s) for s in wf_measure['grim']])
            txt += '\n'
            txt += 'CSM ' + '  '.join(['{:7.3f}'.format(s) for s in wf_measure['csm']])
            txt += '\n'

        output.write(txt)

    def print_wnfsym_irreducible_repr(self, group, axis=None, axis2=None, center=None, output=sys.stdout):

        txt = ''
        for molecule in self._molecules:
            ir_mo = molecule.get_mo_irreducible_representations(group, axis=axis, axis2=axis2, center=center)
            data_wf = molecule.get_wf_irreducible_representations(group, axis=axis, axis2=axis2, center=center)

            # print(data_mo, data_wf)
            txt = '\nMolecule : {}\n'.format(molecule.name)

            sep_line = '     ' + '---------' * len(ir_mo['labels']) + '\n'

            txt += '\nAlpha MOs: Irred. Rep. Decomposition\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_mo['labels']])
            txt += '\n'
            txt += sep_line
            for i, line in enumerate(ir_mo['alpha']):
                txt += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
                txt += '\n'

            txt += '\nBeta MOs: Irred. Rep. Decomposition\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_mo['labels']])
            txt += '\n'
            txt += sep_line
            for i, line in enumerate(ir_mo['beta']):
                txt += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
                txt += '\n'

            txt += '\nWaveFunction: Irred. Rep. Decomposition\n'
            txt += sep_line
            txt += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_mo['labels']])
            txt += '\n'
            txt += sep_line
            txt += 'a-wf' + '  '.join(['{:7.3f}'.format(s) for s in data_wf['alpha']])
            txt += '\n'
            txt += 'b-wf' + '  '.join(['{:7.3f}'.format(s) for s in data_wf['beta']])
            txt += '\n'
            txt += 'WFN ' + '  '.join(['{:7.3f}'.format(s) for s in data_wf['total']])
            txt += '\n'

        output.write(txt)

    def plot_mo_diagram(self, group, axis=None, axis2=None, center=None):

        for molecule in self._molecules:

            ir_mo = molecule.get_mo_irreducible_representations(group, axis=axis, axis2=axis2, center=center)

            ird_a_max = [np.argmax(ird_a_orb) for ird_a_orb in ir_mo['alpha']]
            energies = molecule.electronic_structure.energies

            plt.figure()
            plt.title('{}'.format(molecule.name))

            ax1 = plt.axes()
            ax1.axes.get_xaxis().set_visible(False)  # Hide x axis
            # ax1.axes.get_yaxis().set_visible(True)

            degeneracy = [[energies[0]]]
            for energy in energies[1:]:
                if abs(energy - degeneracy[-1][-1]) < 1e-3:
                    degeneracy[-1].append(energy)
                else:
                    degeneracy.append([energy])

            max_value = 5e-3
            x_center = []
            for ix in degeneracy:
                if len(ix) == 1:
                    x_center.append([0])
                else:
                    x_center.append(np.linspace(-max_value, max_value, len(ix)))
            x_center = [y for x in x_center for y in x]

            plt.scatter(x_center, energies, s=500, marker="_", linewidth=3)
            for i in range(len(energies)):
                plt.text(-max_value * 2, energies[i], ir_mo['labels'][ird_a_max[i]])

        plt.show()

    def plot_sym_energy_evolution(self, group, axis=None, axis2=None, center=None):
        from cosymlib.utils import swap_vectors

        energies = []
        ird_a_max = []
        for idm, molecule in enumerate(self._molecules):

            ir_mo = molecule.get_mo_irreducible_representations(group, axis=axis, axis2=axis2, center=center)

            labels = ir_mo['labels']

            ird_a_max.append(np.array([np.argmax(ird_a_orb) for ird_a_orb in ir_mo['alpha']]))
            energies.append(molecule.electronic_structure.energies)

        energies_x_orbital = np.array(energies).T
        ird_a_x_orbital = np.array(ird_a_max).T

        for i in range(len(ird_a_x_orbital)):
            for j in range(len(ird_a_x_orbital[i])):
                if j == 0:
                    old_ird = ird_a_x_orbital[i][0]
                else:
                    if old_ird != ird_a_x_orbital[i][j]:
                        for k in range(len(ird_a_x_orbital) - i):
                            if old_ird == ird_a_x_orbital[k + i][j]:
                                ird_a_x_orbital[i], ird_a_x_orbital[k + i] = swap_vectors(ird_a_x_orbital[i],
                                                                                          ird_a_x_orbital[k + i], j)
                                energies_x_orbital[i], energies_x_orbital[k + i] = swap_vectors(energies_x_orbital[i],
                                                                                                energies_x_orbital[
                                                                                                    k + i],
                                                                                                j)
                                break
                old_ird = ird_a_x_orbital[i][j]

        for ide, energy in enumerate(energies_x_orbital):
            x = np.arange(len(energy))
            plt.plot(x, energy, marker='_')
            for i in range(len(energy)):
                plt.text(x[i], energy[i] + abs(energy[i]) * 0.001, labels[ird_a_x_orbital[ide][i]])

        plt.show()

    def get_shape_measure(self, label, kind, central_atom=0, fix_permutation=False):
        get_measure = 'get_shape_' + kind
        return [getattr(molecule.geometry, get_measure)(label, central_atom=central_atom,
                                                        fix_permutation=fix_permutation)
                for molecule in self._molecules]

    def get_molecule_path_deviation(self, shape_label1, shape_label2, central_atom=0):
        return [molecule.geometry.get_path_deviation(shape_label1, shape_label2, central_atom) for molecule
                in self._molecules]

    def get_molecule_generalized_coord(self, shape_label1, shape_label2, central_atom=0):
        return [molecule.geometry.get_generalized_coordinate(shape_label1, shape_label2, central_atom)
                for molecule in self._molecules]

    # TODO: This may be placed inside Shape class
    def get_path_parameters(self, shape_label1, shape_label2, central_atom=0):

        if type(shape_label1) is Geometry:
            label1 = shape_label1
            label1_name = shape_label1.name
        else:
            label1 = shape_label1
            label1_name = shape_label1
        if type(shape_label2) is Geometry:
            label2 = shape_label2
            label2_name = shape_label2.name
        else:
            label2 = shape_label2
            label2_name = shape_label2

        csm = {label1_name: self.get_shape_measure(label1, 'measure', central_atom),
               label2_name: self.get_shape_measure(label2, 'measure', central_atom)}
        devpath = self.get_molecule_path_deviation(label1, label2, central_atom)
        generalized_coord = self.get_molecule_generalized_coord(label1, label2, central_atom)

        return csm, devpath, generalized_coord

    def print_minimum_distortion_path_shape(self, shape_label1, shape_label2, central_atom=0,
                                            min_dev=0, max_dev=15, min_gco=0, max_gco=101,
                                            num_points=20, output_name=None):

        if output_name is not None:
            output = open(output_name + '_pth.csv', 'w')
            output2 = open(output_name + '_pth.xyz', 'w')
            output3 = open(output_name, 'w')
        else:
            output = sys.stdout
            output2 = sys.stdout
            output3 = sys.stdout

        csm, devpath, gen_coord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom)

        txt_params = 'Deviation threshold to calculate Path deviation function: {:2.1f}% - {:2.1f}%\n'.format(min_dev, max_dev)
        txt_params += 'Deviation threshold to calculate Generalized Coordinate: {:2.1f}% - {:2.1f}%\n'.format(min_gco, max_gco)
        txt_params += '\n'
        txt_params += '{:9} '.format('structure'.upper())
        for csm_label in list(csm.keys()):
            txt_params += '{:^8} '.format(csm_label)
        txt_params += '{:^8} {:^8}'.format('DevPath', 'GenCoord')
        txt_params += '\n'

        filter_mask = [min_dev < dv < max_dev and min_gco < gc < max_gco for dv, gc in zip(devpath, gen_coord)]

        for idx, molecule in enumerate(self._molecules):
            if not filter_mask[idx]:
                continue

            txt_params += '{:9} '.format(molecule.name.strip())
            for label in list(csm.keys()):
                txt_params += '{:^8.3f} '.format(csm[label][idx])
            txt_params += '{:^8.1f} {:^8.1f}'.format(devpath[idx], gen_coord[idx])
            txt_params += '\n'

        txt_params += 'skipped {} structure/s\n\n'.format(filter_mask.count(False))
        output3.write(txt_params)

        if isinstance(shape_label1, Geometry):
            label1_name = shape_label1.name
        else:
            label1_name = shape_label1
        if isinstance(shape_label2, Geometry):
            label2_name = shape_label2.name
        else:
            label2_name = shape_label2

        path = get_shape_path(shape_label1, shape_label2, num_points)
        txt_path = 'Minimum distortion path\n'
        txt_path += ' {:^6}  {:^6}\n'.format(label1_name, label2_name)
        for idx, value in enumerate(path[0]):
            txt_path += '{:6.3f}  {:6.3f}'.format(path[0][idx], path[1][idx])
            txt_path += '\n'
        txt_path += '\n'
        output.write(txt_path)

        test_structures = []
        for ids, structure in enumerate(path[2]):
            test_structures.append(Geometry(symbols=['' for _ in range(len(structure))],
                                            positions=structure, name='map_structure{}'.format(ids)))
        output2.write(file_io.get_file_xyz_txt(test_structures))

        if output_name is None:
            import matplotlib.pyplot as plt
            plt.plot(path[0], path[1], 'k', linewidth=2.0)
            plt.scatter(np.array(csm[label1_name])[filter_mask],
                        np.array(csm[label2_name])[filter_mask], linewidths=0.01)
            plt.xlabel(label1_name)
            plt.ylabel(label2_name)
            plt.show()

    def get_point_group(self, tol=0.01):
        return [molecule.geometry.get_pointgroup(tol) for molecule in self._molecules]
