from symeess import file_io
from itertools import compress
__version__ = 0.61


class Symeess:

    def __init__(self, file_name, shape_label=None, shape_options=None, output_name=None):

        self._file_name = file_name
        self._shape_label = shape_label.split()
        self._results = {}
        self._molecules = self.read_input()

        # Shape measures
        if shape_options is not None:
            self.shape(shape_options, output_name)

    def read_input(self):
        return file_io.read(self._file_name)

    def shape(self, shape_options, output_name):

        shape_type_measures = ['get_ideal_structure', 'get_measure', 'get_test_structure']
        shape_choices = list(compress(shape_type_measures, shape_options[:-1]))

        if shape_choices:
            for key, molecule in self._molecules.items():
                self._results[key] = {}
                self._results[key]['symbols'] = molecule.geometry.get_symbols()
                for label in self._shape_label:
                    self._results[key][label] = {}
                    for pattern in shape_choices:
                        name = pattern.replace('get_', '')
                        self._results[key][label][name] = (getattr(molecule.geometry, pattern)
                                                        (shape_label=label, central_atom=shape_options[-1]))

            # Canviar el metode d'escriptura
            file_io.write_shape(output_name, self._results, self._shape_label, shape_choices)