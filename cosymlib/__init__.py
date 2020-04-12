__version__ = '0.7.4'

from cosymlib.molecule import Molecule, Geometry
from cosymlib import file_io
from cosymlib import tools
from cosymlib.file_io.shape import get_shape_measure_data_txt, write_minimal_distortion_path_analysis, write_shape_map
from cosymlib.utils import get_shape_map, plot_molecular_orbital_diagram, plot_symmetry_energy_evolution
from cosymlib.shape.tools import get_structure_references

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
                 ignore_connectivity=False):

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
            if ignore_connectivity:
                molecule.geometry.set_connectivity(None)

    def get_n_atoms(self):
        n_atoms_unique_list = np.unique([mol.geometry.get_n_atoms() for mol in self._molecules])
        if len(n_atoms_unique_list) > 1:
            raise Exception('Not all structures have same number of atoms')

        return n_atoms_unique_list[0]

    def get_geometries(self):
        return [mol.geometry for mol in self._molecules]

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
        shape_results_measures = []
        references_names = []
        for reference in reference_list:
            shape_results_measures.append(self.get_shape_measure(reference, 'measure', central_atom,
                                                                 fix_permutation))
            if type(reference) is Geometry:
                references_names.append(reference.name)
            else:
                references_names.append(reference)

        output.write(get_shape_measure_data_txt(shape_results_measures, molecules_names, references_names))

    def print_shape_structure(self, shape_reference, central_atom=0, fix_permutation=False, output=sys.stdout):
        """
        Method that prints to file shape's structure
        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                                Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """

        shape_results_structures = []
        references = []
        for reference in shape_reference:
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
                geometries.append(Geometry(symbols=molecule.geometry.get_symbols(),
                                           positions=shape_results_structures[idl][idm],
                                           name=molecule.name + '_' + reference))

        for geometry in geometries:
            output.write(file_io.get_file_xyz_txt(geometry))

    def print_path_parameters(self, shape_label1, shape_label2, central_atom=0,
                              maxdev=15, mindev=0, maxgco=101, mingco=0, output=sys.stdout):

        #if output_name is not None:
        #    output = open(output_name + '_tab.csv', 'w')
        #else:
        #    output = sys.stdout

        csm, devpath, GenCoord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom,
                                                          maxdev=maxdev, mindev=mindev, maxgco=maxgco, mingco=mingco)
        names_order = [molecule.name for molecule in self._molecules]
        txt = write_minimal_distortion_path_analysis(csm, devpath, GenCoord, maxdev, mindev,
                                                     mingco, maxgco, names_order)
        output.write(txt)

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

            txt += 'Symmetry measure {:.5f}\n'.format(molecule.geometry.get_symmetry_measure(kwargs))
            txt += sep_line

        output.write(txt)

    def print_geometric_symmetry_measure(self, label, multi=1, central_atom=0, connect_thresh=1.1, center=None, output=sys.stdout):
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

    def print_symmetry_nearest_structure(self, label, multi=1, central_atom=0, connect_thresh=1.1, center=None, output=sys.stdout):
        kwargs = _get_symgroup_arguments(locals())

        for idm, molecule in enumerate(self._molecules):
            geometry = Geometry(symbols=molecule.geometry.get_symbols(),
                                positions=molecule.geometry.get_symmetry_nearest_structure(**kwargs),
                                 name=molecule.name)

            output.write(file_io.get_file_xyz_txt(geometry))

    def print_wnfsym_measure_verbose(self, group,
                                     vector_axis1=None,
                                     vector_axis2=None,
                                     center=None,
                                     output=sys.stdout,
                                     n_molecule=0):

        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        txt = file_io.symmetry.get_operated_matrices_txt(group, self._molecules[n_molecule],
                                                               wfnsym_results[n_molecule])
        txt += file_io.symmetry.get_overlap_analysis_txt(wfnsym_results[n_molecule])
        txt += file_io.symmetry.get_ir_analysis_txt(wfnsym_results[n_molecule])
        output.write(txt)

    def print_wnfsym_sym_matrices(self, group, vector_axis1=None, vector_axis2=None,
                                  center=None, n_molecule=0, output=sys.stdout):

        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        txt = file_io.symmetry.get_operated_matrices_txt(group, self._molecules[0],
                                                               wfnsym_results[n_molecule])
        output.write(txt)

    def print_wnfsym_sym_ovelap(self, group, vector_axis1=None, vector_axis2=None,
                                center=None, n_molecule=0, output=sys.stdout):

        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        txt = file_io.symmetry.get_overlap_analysis_txt(wfnsym_results[n_molecule])
        output.write(txt)

    def print_wnfsym_ireducible_repr(self, group, vector_axis1=None, vector_axis2=None, center=None,
                                     n_molecule=0, output=sys.stdout):

        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        txt = file_io.symmetry.get_ir_analysis_txt(wfnsym_results[n_molecule])
        output.write(txt)

    def plot_mo_diagram(self, group, vector_axis1=None, vector_axis2=None, center=None, n_molecule=0):
        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        plot_molecular_orbital_diagram(self._molecules[0], wfnsym_results[n_molecule])

    def plot_sym_energy_evolution(self, group, vector_axis1=None, vector_axis2=None, center=None):
        wfnsym_results = self._get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        plot_symmetry_energy_evolution(self._molecules, wfnsym_results)

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

    def get_path_parameters(self, shape_label1, shape_label2, central_atom=0, maxdev=15, mindev=0,
                            maxgco=101, mingco=0):

        def filter_results(results, criteria, max, min):
            filter = [True if x <= max and x >= min else False for x in criteria]
            results = [i for indx, i in enumerate(results) if filter[indx] == True]
            return results

        if type(shape_label1) is Geometry:
            label1 = shape_label1.get_positions()
            label1_name = shape_label1.name
        else:
            label1 = shape_label1
            label1_name = shape_label1
        if type(shape_label2) is Geometry:
            label2 = shape_label2.get_positions()
            label2_name = shape_label2.name
        else:
            label2 = shape_label2
            label2_name = shape_label2

        csm = {label1_name: self.get_shape_measure(label1, 'measure', central_atom),
               label2_name: self.get_shape_measure(label2, 'measure', central_atom)}
        devpath = self.get_molecule_path_deviation(label1, label2, central_atom)
        generalized_coord = self.get_molecule_generalized_coord(label1, label2, central_atom)
        criteria = devpath
        devpath = filter_results(devpath, criteria, maxdev, mindev)
        generalized_coord = filter_results(generalized_coord, criteria, maxdev, mindev)
        criteria = generalized_coord
        devpath = filter_results(devpath, criteria, maxgco, mingco)
        generalized_coord = filter_results(generalized_coord, criteria, maxgco, mingco)
        return csm, devpath, generalized_coord

    def _get_wfnsym_results(self, group, vector_axis1, vector_axis2, center):
        return [molecule.get_mo_symmetry(group, vector_axis1=vector_axis1, vector_axis2=vector_axis2, center=center)
                for molecule in self._molecules]

    def print_minimum_distortion_path_shape(self, shape_label1, shape_label2, central_atom=0,
                                            num_points=20, output_name=None):

        if output_name is not None:
            output = open(output_name + '_pth.csv', 'w')
            output2 = open(output_name + '_pth.xyz', 'w')
            output3 = open(output_name, 'w')
        else:
            output = sys.stdout
            output2 = sys.stdout
            output3 = sys.stdout

        self.print_path_parameters(shape_label1, shape_label2, central_atom=central_atom, output=output3)

        if type(shape_label1) is Geometry:
            label1 = shape_label1.get_positions()
            label1_name = shape_label1.name
        else:
            label1 = shape_label1
            label1_name = shape_label1
        if type(shape_label2) is Geometry:
            label2 = shape_label2.get_positions()
            label2_name = shape_label2.name
        else:
            label2 = shape_label2
            label2_name = shape_label2

        path_parameters = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom)[0]
        path = get_shape_map(label1, label2, num_points)
        output.write(write_shape_map(label1_name, label2_name, path))
        test_structures = []
        for ids, structure in enumerate(path[2]):
            test_structures.append(Geometry(symbols=['' for _ in range(len(structure))],
                                            positions=structure, name='map_structure{}'.format(ids)))
        output2.write(file_io.get_file_xyz_txt(test_structures))
        if output_name is None:
            import matplotlib.pyplot as plt
            plt.plot(path[0], path[1], 'k', linewidth=2.0)
            plt.scatter(path_parameters[label1_name], path_parameters[label2_name], linewidths=0.01)
            plt.xlabel(label1_name)
            plt.ylabel(label2_name)
            plt.show()

    def get_point_group(self, tol=0.01):
        return [molecule.geometry.get_pointgroup(tol) for molecule in self._molecules]
