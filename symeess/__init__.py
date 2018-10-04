from symeess.molecule import Molecule
from symeess import file_io
from symeess.file_io import shape2file
from symeess.shape import maps
import matplotlib.pyplot as plt
__version__ = 0.6


class Symeess:
    """
    Main class of symeess program that perform all the jobs
    """

    def __init__(self):

        self._molecules = []

    def set_molecules(self, molecules):
        if isinstance(molecules, list):
            for molecule in molecules:
                try:
                    molecule.geometry.get_positions()
                    self._molecules.append(molecule)
                except AttributeError:
                    try:
                        molecule.get_positions()
                        self._molecules.append(Molecule(molecule))
                    except AttributeError:
                        raise AttributeError('Molecule object not found')
        else:
            self._molecules.append(molecules)

    def read_molecules(self, input_file):
        geometries = file_io.read_input_file(input_file)
        self.set_molecules(geometries)

    def write_shape_measure_2file(self, shape_label, central_atom=None, output_name='symeess_shape'):
        """
        Method that prints to file shape's measure
        :param shape_label: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's measure in the output_name.tab file
        """
        shape_results_measures = [self.get_shape_measure(label, 'measure', central_atom) for label in shape_label]
        molecules_name = [molecule.get_name() for molecule in self._molecules]
        shape2file.write_shape_measure_data(shape_results_measures, molecules_name, shape_label, output_name=output_name)

    def write_shape_structure_2file(self, shape_label, central_atom=None, output_name='symeess_shape'):
        """
        Method that prints to file shape's structure
        :param shape_label: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """
        initial_geometry = [molecule.geometry.get_positions() for molecule in self._molecules]
        shape_results_structures = [self.get_shape_measure(label, 'structure', central_atom) for label in shape_label]
        molecules_name = [molecule.get_name() for molecule in self._molecules]
        symbols = [molecule.geometry.get_symbols() for molecule in self._molecules]
        shape_results_measures = [self.get_shape_measure(label, 'measure', central_atom) for label in shape_label]
        shape2file.write_shape_structure_data(initial_geometry, shape_results_structures, shape_results_measures,
                                              symbols, molecules_name, shape_label,
                                              output_name=output_name)

    def write_path_parameters_2file(self, shape_label1, shape_label2, central_atom=None,
                                    maxdev=15, mindev=0, maxgco=101, mingco=0, output_name='symeess_shape'):

        shape, devpath, GenCoord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom,
                                                            maxdev=maxdev, mindev=mindev, maxgco=maxgco, mingco=mingco)
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape2file.write_minimal_distortion_path_analysis(shape_label1, shape_label2, shape, devpath, GenCoord,
                                                          maxdev, mindev, mingco, maxgco, names_order,
                                                          output_name='symeess_shape')

    def write_symgroup_measure(self, group, multi=1, central_atom=None, output_name='symeess_symgroup'):
        results = self.get_symgroup_measure(group=group, multi=multi, central_atom=central_atom)
        file_io.write_symgroup_measure(group, [molecule.geometry for molecule in self._molecules], results, output_name)

    def write_wnfsym_measure_2file(self, label, vector_axis1, vector_axis2, center, output_name='symeess_wfnsym'):
        wfnsym_results = self.get_wfnsym_measure(label, vector_axis1, vector_axis2, center)
        file_io.write_wfnsym_measure(label, self._molecules[0].geometry, wfnsym_results, output_name)

    def get_shape_measure(self, label, type, central_atom=None):
        """
        :param label: reference polyhedra label which user will compare with his polyhedra.
                      Reference labels can be found in [#f1]_
        :param type: type of measure that the user is going to use
        :param central_atom: position of the central atom in molecule if exist
        :return:
        """
        get_measure = 'get_shape_' + type
        shape_data = [getattr(molecule.geometry, get_measure)(label, central_atom=central_atom)
                 for molecule in self._molecules]
        return shape_data

    def get_molecule_path_deviation(self, shape_label1, shape_label2, central_atom=None):
        path_deviation = [molecule.geometry.get_path_deviation(shape_label1, shape_label2, central_atom) for molecule
                          in self._molecules]
        return path_deviation

    def get_molecule_GenCoord(self, shape_label1, shape_label2, central_atom=None):
        GenCoord = [molecule.geometry.get_generalized_coordinate(shape_label1, shape_label2, central_atom)
                    for molecule in self._molecules]
        return GenCoord

    def get_path_parameters(self, shape_label1, shape_label2, central_atom=None, maxdev=15, mindev=0,
                            maxgco=101, mingco=0):

        csm = {shape_label1: self.get_shape_measure(shape_label1, 'measure', central_atom),
               shape_label2: self.get_shape_measure(shape_label2, 'measure', central_atom)}
        devpath = self.get_molecule_path_deviation(shape_label1, shape_label2, central_atom)
        generalized_coord = self.get_molecule_GenCoord(shape_label1, shape_label2, central_atom)
        criteria = devpath
        devpath = self.get_filtered_results(devpath, criteria, maxdev, mindev)
        generalized_coord = self.get_filtered_results(generalized_coord, criteria, maxdev, mindev)
        criteria = generalized_coord
        devpath = self.get_filtered_results(devpath, criteria, maxgco, mingco)
        generalized_coord = self.get_filtered_results(generalized_coord, criteria, maxgco, mingco)
        return csm, devpath, generalized_coord

    def get_filtered_results(self, results, criteria, max, min):
        pathdev_filter = [True if x <= max and x >= min else False for x in criteria]
        results = [i for indx, i in enumerate(results) if pathdev_filter[indx] == True]
        return results

    def get_symgroup_measure(self, group, multi=1, central_atom=None):
        results = [molecule.geometry.get_symmetry_measure(label=group, multi=multi, central_atom=central_atom) for
                   molecule in self._molecules]
        return results

    def get_wfnsym_measure(self, label, vector_axis1, vector_axis2, center):
        results = self._molecules[0].get_mo_symmetry(label,
                                                     vector_axis2=vector_axis2,
                                                     vector_axis1=vector_axis1,
                                                     center=center)
        return results


def write_minimum_distortion_path_shape_2file(shape_label1, shape_label2, central_atom=None, num_points=50, show=False):
    path = get_shape_map(shape_label1, shape_label2, central_atom, num_points)
    shape2file.write_shape_map(shape_label1, shape_label2, path)
    if show:
        plt.plot(path[0], path[1], linewidth=2.0)
        plt.xlabel(shape_label1)
        plt.ylabel(shape_label2)
        plt.show()


def get_shape_map(shape_label1, shape_label2, central_atom, num_points):
    x, y = maps.get_shape_map(shape_label1, shape_label2, central_atom, num_points)
    return x, y