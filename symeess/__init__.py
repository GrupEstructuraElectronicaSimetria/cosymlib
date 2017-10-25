from symeess import file_io
__version__ = 0.7


class Symeess:

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

    def write_shape_measure_2file(self, shape_label, central_atom=None, output_name='../examples/symeess_shape'):
        """

        :param shape_label:
        :param central_atom: position of central atom if exist (default None)
        :param output_name: custom name without extension (default ../examples/symeess_shape)
        :return:
        """
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape = {}
        for label in shape_label.split():
            shape[label] = [molecule.geometry.get_shape_measure(label, central_atom=central_atom)
                            for molecule in self._molecules]
        file_io.write_shape_data(shape, shape_label.split(), names_order, 'measure', output_name)

    def write_shape_structure_2file(self, shape_label, central_atom=None, output_name='../examples/symeess_shape'):
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape = {}
        for label in shape_label.split():
            shape[label] = [molecule.geometry.get_shape_structure(label, central_atom=central_atom)
                            for molecule in self._molecules]
        shape['symbols'] = [molecule.geometry.get_symbols() for molecule in self._molecules]
        file_io.write_shape_data(shape, shape_label.split(), names_order, 'structure', output_name)
