from symeess.shape import measure, structure_measure, test_structure


class Geometry:
    def __init__(self, structure=None):

        self._symbols = []
        self._positions = []
        self._c_positions = []

        for elements in structure:
            try:
                int(elements[0][1])
                self._symbols.append(elements[0][0])
            except (ValueError, TypeError, IndexError):
                self._symbols.append(elements[0][:2])
            self._positions.append([float(j) for j in elements[1:]])

        self._n_atoms = None
        self._shape_ideal = 0
        self._central_atom = None
        self._shape_test_structure = []
        self._shape_measures = {}

    def get_positions(self):
        return self._positions

    def get_n_atoms(self):
        self._n_atoms = len(self._positions)
        return self._n_atoms

    def get_symbols(self):
        return self._symbols

    def set_measure(self, shape_label, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_label
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_measures[shape_label]['measure'] = measure(self.get_positions(), n_atoms,
                                                               self._shape_ideal, self._central_atom)

    def set_ideal_structure(self, shape_label, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_label
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_measures[shape_label]['measure'], self._shape_measures[shape_label]['structure'] = \
            structure_measure(self.get_positions(), n_atoms, self._shape_ideal, self._central_atom)

    def set_test_structure(self, shape_label, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_label
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_test_structure = test_structure(self.get_positions(), n_atoms,
                                                    self._shape_ideal, self._central_atom)

    def get_measure(self, shape_label, central_atom=None):
        if shape_label not in self._shape_measures:
            self._shape_measures[shape_label] = {}
        if 'measure' not in self._shape_measures[shape_label]:
            self.set_measure(shape_label, central_atom)
        return self._shape_measures[shape_label]['measure']

    def get_ideal_structure(self, shape_label, central_atom=None):
        if shape_label not in self._shape_measures:
            self._shape_measures[shape_label] = {}
        if 'structure' not in self._shape_measures[shape_label]:
            self.set_ideal_structure(shape_label, central_atom)
        return self._shape_measures[shape_label]['structure']

    def get_test_structure(self, shape_label, central_atom=None):
        self.set_test_structure(shape_label, central_atom)
        return self._shape_test_structure
