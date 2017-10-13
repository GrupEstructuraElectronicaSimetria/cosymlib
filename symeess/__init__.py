from symeess import file_io
from itertools import compress
__version__ = 0.1


class Symeess:

    def __init__(self, file_name, shape_label=None, shape_options=None, output_name=None):

        self._file_name = file_name
        self._shape_label = shape_label
        self._results = {}
        self._molecule = self.read_input()

        # Shape measures
        shape_type_measures = ['get_ideal_structure', 'get_measure', 'get_test_structure']
        shape_choices = list(compress(shape_type_measures, shape_options[:-1]))
        if shape_choices:
            for pattern in shape_choices:
                self._results[pattern] = (getattr(self._molecule.geometry, pattern)(shape_reference=self._shape_label,
                                                                                    central_atom=shape_options[-1]))
            file_io.write(output_name, self._results, self._shape_label, shape_choices)

    def read_input(self):
        return file_io.read(self._file_name)
