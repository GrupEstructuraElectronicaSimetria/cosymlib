from symeess import file_io
__version__ = 0.1


class Symeess:
    """
    Main class of symeess program that perform all the jobs

    """

    def __init__(self, input_data):

        self._shape_label = None
        self._results = []
        self._molecules = input_data

    # def set_shape_tstructure(self, shape_label, central_atom):
    #     if self._shape_label is None:
    #         self._shape_label = shape_label.split()
    #     for idx, molecule in enumerate(self._molecules):
    #         if len(self._results) == idx:
    #             self._results.append({})
    #         for label in self._shape_label:
    #             if label not in self._results[idx]:
    #                 self._results[idx][label] = {}
    #             self._results[idx][label]['test_structure'] = molecule.geometry.get_test_structure(label, central_atom)

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
            shape[label] = [molecule.geometry.get_shape_measure(label, central_atom=central_atom)
                            for molecule in self._molecules]
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
            shape[label] = [molecule.geometry.get_shape_structure(label, central_atom=central_atom)
                            for molecule in self._molecules]
        shape['symbols'] = [molecule.geometry.get_symbols() for molecule in self._molecules]
        file_io.write_shape_data(shape, shape_label.split(), names_order, 'structure', output_name)
