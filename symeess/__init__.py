from symeess import file_io
__version__ = 0.61


class Symeess:

    def __init__(self, file_name):

        self._file_name = file_name
        self._shape_label = None
        self._results = {}
        self._molecules = self.read_input()

    def read_input(self):
        return file_io.read(self._file_name)

    def set_shape_measure(self, shape_label, central_atom):
        if self._shape_label is None:
            self._shape_label = shape_label.split()
        for key, molecule in self._molecules.items():
            if not self._results[key]:
                self._results[key] = {}
            for label in self._shape_label:
                if not self._results[key][label]:
                    self._results[key][label] = {}
                self._results[key][label]['measure'] = molecule.geometry.get_measure(label, central_atom)

    def set_shape_structure(self, shape_label, central_atom):
        if self._shape_label is None:
            self._shape_label = shape_label.split()
        for key, molecule in self._molecules.items():
            self._results[key] = {}
            self._results[key]['symbols'] = molecule.geometry.get_symbols()
            for label in self._shape_label:
                self._results[key][label] = {}
                self._results[key][label]['ideal_structure'] = molecule.geometry.get_ideal_structure(label,
                                                                                                     central_atom)

    def set_shape_tstructure(self, shape_label, central_atom):
        if self._shape_label is None:
            self._shape_label = shape_label.split()
        for key, molecule in self._molecules.items():
            if not self._results[key]:
                self._results[key] = {}
                self._results[key]['symbols'] = molecule.geometry.get_symbols()
            for label in self._shape_label:
                if not self._results[key][label]:
                    self._results[key][label] = {}
                self._results[key][label]['test_structure'] = molecule.geometry.get_test_structure(label, central_atom)

    def write_shape(self, output_name):
        file_io.write_shape(output_name, self._results, self._shape_label)
