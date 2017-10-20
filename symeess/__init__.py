from symeess import file_io
__version__ = 0.61


class Symeess:

    def __init__(self, input_data):

        self._shape_label = None
        self._results = []
        self._molecules = input_data

    def _calculate_shape(self, shape_label, central_atom, measure='measure'):
        get_measure = 'get_shape_'+measure
        if self._shape_label is None:
            self._shape_label = shape_label.split()
        for idx, molecule in enumerate(self._molecules):
            if len(self._results) == idx:
                self._results.append({})
                self._results[idx]['symbols'] = molecule.geometry.get_symbols()
            for label in self._shape_label:
                if label not in self._results[idx]:
                    self._results[idx][label] = {}
                self._results[idx][label][measure] = getattr(molecule.geometry, get_measure)(label, central_atom)

    def set_shape_tstructure(self, shape_label, central_atom):
        if self._shape_label is None:
            self._shape_label = shape_label.split()
        for idx, molecule in enumerate(self._molecules):
            if len(self._results) == idx:
                self._results.append({})
            for label in self._shape_label:
                if label not in self._results[idx]:
                    self._results[idx][label] = {}
                self._results[idx][label]['test_structure'] = molecule.geometry.get_test_structure(label, central_atom)

    def write_shape_measure(self, shape_label, central_atom=None, output_name='../examples/symeess_shape'):
        self._calculate_shape(shape_label, central_atom, measure='measure')
        names_order = [molecule.get_name() for molecule in self._molecules]
        file_io.write_shape_data(output_name, self._results, self._shape_label, names_order, 'measure')

    def write_shape_structure(self, shape_label, central_atom=None, output_name='../examples/symeess_shape'):
        self._calculate_shape(shape_label, central_atom, measure='structure')
        names_order = [molecule.get_name() for molecule in self._molecules]
        file_io.write_shape_data(output_name, self._results, self._shape_label, names_order, 'structure')
