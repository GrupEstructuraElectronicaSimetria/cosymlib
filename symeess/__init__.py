import symeess.file_io as file_io
from symeess import maps
import matplotlib.pyplot as plt
__version__ = 0.1


class Symeess:
    """
    Main class of symeess program that perform all the jobs

    """

    def __init__(self, input_data):

        self._molecules = input_data

    def get_shape_measure(self, label, type, central_atom=None):
        """

        :param label: reference polyhedra label which user will compare with his polyhedra.
                      Reference labels can be found in [#f1]_
        :param type: type of measure that the user is going to use
        :param central_atom: position of the central atom in molecule if exist
        :return:
        """
        get_measure = 'get_shape_' + type
        shape = [getattr(molecule.geometry, get_measure)(label, central_atom=central_atom)
                 for molecule in self._molecules]
        return shape

    def write_shape_measure_2file(self, shape_label, central_atom=None, output_name='symeess_shape'):
        """
        Method that prints to file shape's measure

        :param shape_label: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's measure in the output_name.tab file
        """
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape = {}
        for label in shape_label.split():
            shape[label] = self.get_shape_measure(label, 'measure', central_atom)
        file_io.write_shape_data(shape, shape_label.split(), names_order, 'measure', output_name)

    def write_shape_structure_2file(self, shape_label, central_atom=None, output_name='symeess_shape'):
        """
        Method that prints to file shape's structure

        :param shape_label: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape = {}
        for label in shape_label.split():
            shape[label] = self.get_shape_measure(label, 'structure', central_atom)
        shape['symbols'] = [molecule.geometry.get_symbols() for molecule in self._molecules]
        file_io.write_shape_data(shape, shape_label.split(), names_order, 'structure', output_name)

    def get_molecule_path_deviation(self, shape_label1, shape_label2, central_atom=None):
        path_deviation = [molecule.geometry.get_path_deviation(shape_label1, shape_label2, central_atom) for molecule
                          in self._molecules]
        return path_deviation

    def minimum_distortion_path_shape_2file(self, shape_label1, shape_label2, central_atom=None, num_points=50, show=False):
        path = self.get_shape_map(shape_label1, shape_label2, central_atom, num_points)
        file_io.write_shape_map_2file(shape_label1, shape_label2, path)
        if show:
            plt.plot(path[0], path[1], linewidth=2.0)
            plt.xlabel(shape_label1)
            plt.ylabel(shape_label2)
            plt.show()


    def get_shape_map(self, shape_label1, shape_label2, central_atom, num_points):
        x, y = maps.get_shape_map(shape_label1, shape_label2, central_atom, num_points)
        return x, y
