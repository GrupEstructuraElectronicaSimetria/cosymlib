from cosymlib import shape, tools
from cosymlib.symmetry import symgroup
import numpy as np


class Geometry:
    def __init__(self,
                 symbols=[],
                 positions=None,
                 name=None):

        self._central_atom = None
        self._symbols = []
        self._positions = []
        self._atom_groups = list(symbols)

        if name is not None:
            self._name = name
        else:
            self._name = ' ' * 5

        for symbol in symbols:
            try:
                int(symbol)
                self._symbols.append(tools.atomic_number_to_element(int(symbol)))
            except (ValueError, TypeError):
                self._symbols.append(symbol.capitalize())
                for ida, a in enumerate(symbol):
                    try:
                        int(a)
                        self._symbols[-1] = self._symbols[-1][:ida]
                        break
                    except (ValueError, TypeError, IndexError):
                        pass

        try:
            float(positions[0])
            for symbol in positions:
                self._positions.append(float(symbol))
            self._positions = list(chunks(self._positions, 3))
        except (ValueError, TypeError, IndexError):
            for symbol in positions:
                self._positions.append([float(j) for j in symbol])

        self._positions = np.array(self._positions)
        self._shape = shape.Shape(self)
        self._symgroup = symgroup.Symgroup(self)

    def get_name(self):
        return self._name

    def set_name(self, name):
        self._name = name

    def get_positions(self):
        return self._positions

    def get_n_atoms(self):
        return len(self.get_positions())

    def get_symbols(self):
        return self._symbols

    def get_shape_measure(self, shape_label, central_atom=0, fix_permutation=False):
        return self._shape.measure(shape_label, central_atom=central_atom, fix_permutation=fix_permutation)

    def get_shape_structure(self, shape_label, central_atom=0, fix_permutation=False):
        return self._shape.structure(shape_label, central_atom=central_atom, fix_permutation=fix_permutation)

    def get_symmetry_measure(self, label, central_atom=0, multi=1, symbols=True):
        return self._symgroup.results(label, central_atom=central_atom, multi=multi)

    def get_path_deviation(self, shape_label1, shape_label2, central_atom=0):
        return self._shape.get_path_deviation(shape_label1, shape_label2, central_atom)

    def get_generalized_coordinate(self, shape_label1, shape_label2, central_atom=0):
        return self._shape.get_generalized_coordinate(shape_label1, shape_label2, central_atom)


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]
