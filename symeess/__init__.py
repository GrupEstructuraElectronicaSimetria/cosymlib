import symeess.file_io as file_io
from symeess.file_io import shape2file
from symeess.shape import maps
import matplotlib.pyplot as plt
__version__ = 0.1


class Symeess:
    """
    Main class of symeess program that perform all the jobs

    """

    def __init__(self, input_data):

        self._molecules = input_data


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
        shape2file.write_shape_data(shape, shape_label.split(), names_order, 'measure', output_name)

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
        shape2file.write_shape_data(shape, shape_label.split(), names_order, 'structure', output_name)

    def write_path_parameters_2file(self, shape_label1, shape_label2, central_atom=None,
                                    maxdev=15, mindev=0, maxgco=101, mingco=0):
        """

        :param shape_label1:
        :param shape_label2:
        :param central_atom:
        :param maxdev:
        :param mindev:
        :param maxgco:
        :param mingco:
        :return:
        """
        shape, devpath, GenCoord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom,
                                                            maxdev=maxdev, mindev=mindev, maxgco=maxgco, mingco=mingco)
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape2file.write_minimal_distortion_path_analysis(shape_label1, shape_label2, shape, devpath, GenCoord,
                                                          maxdev, mindev, mingco, maxgco, names_order,
                                                          output_name='symeess_shape')


    def write_minimum_distortion_path_shape_2file(self, shape_label1, shape_label2, central_atom=None,
                                                  num_points=50, show=False):
        path = self.get_shape_map(shape_label1, shape_label2, central_atom, num_points)
        shape2file.write_shape_map(shape_label1, shape_label2, path)
        if show:
            plt.plot(path[0], path[1], linewidth=2.0)
            plt.xlabel(shape_label1)
            plt.ylabel(shape_label2)
            plt.show()

    def write_wyfsym_measure_2file(self, label, VAxis1, VAxis2, RCread, output_name='symeess_wyfsym'):
        results = self.get_wyfsym_measure(label, VAxis1, VAxis2, RCread)
        file_io.write_wyfsym_measure(results, output_name)

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

    def get_molecule_path_deviation(self, shape_label1, shape_label2, central_atom=None):
        path_deviation = [molecule.geometry.get_path_deviation(shape_label1, shape_label2, central_atom) for molecule
                          in self._molecules]
        return path_deviation

    def get_molecule_GenCoord(self, shape_label1, shape_label2, central_atom=None):
        GenCoord = [molecule.geometry.get_generalized_coordinate(shape_label1, shape_label2, central_atom)
                    for molecule in self._molecules]
        return GenCoord

    def get_shape_map(self, shape_label1, shape_label2, central_atom, num_points):
        x, y = maps.get_shape_map(shape_label1, shape_label2, central_atom, num_points)
        return x, y

    def get_path_parameters(self, shape_label1, shape_label2, central_atom=None, maxdev=15, mindev=0,
                            maxgco=101, mingco=0):
        shape = {}
        shape[shape_label1] = self.get_shape_measure(shape_label1, 'measure', central_atom)
        shape[shape_label2] = self.get_shape_measure(shape_label2, 'measure', central_atom)
        devpath = self.get_molecule_path_deviation(shape_label1, shape_label2, central_atom)
        GenCoord = self.get_molecule_GenCoord(shape_label1, shape_label2, central_atom)
        criteria = devpath
        devpath = self.get_filtered_results(devpath, criteria, maxdev, mindev)
        GenCoord = self.get_filtered_results(GenCoord, criteria, maxdev, mindev)
        criteria = GenCoord
        devpath = self.get_filtered_results(devpath, criteria, maxgco, mingco)
        GenCoord = self.get_filtered_results(GenCoord, criteria, maxgco, mingco)
        return shape, devpath, GenCoord

    def get_filtered_results(self, results, criteria, max, min):
        pathdev_filter = [True if x <= max and x >= min else False for x in criteria]
        results = [i for indx, i in enumerate(results) if pathdev_filter[indx] == True]
        return results

    def get_wyfsym_measure(self, label, VAxis1, VAxis2, RCread):
        results = [molecule.electronic_structure.get_wyfsym_measure(8, 3, VAxis1, VAxis2, RCread)
                   for molecule in self._molecules]
        return results[0]
