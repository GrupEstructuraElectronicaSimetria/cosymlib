from symeess import shape


class Geometry:
    def __init__(self, structure=None):

        self._symbols = []
        self._positions = []

        for elements in structure:
            try:
                int(elements[0][1])
                self._symbols.append(elements[0][0].capitalize())
            except (ValueError, TypeError, IndexError):
                self._symbols.append(elements[0][:2].capitalize())
            self._positions.append([float(j) for j in elements[1:]])

        self._n_atoms = None
        self._shape_label = 0
        self._central_atom = None
        self._shape_measures = {}
        self._path_deviation = {}
        self._GenCoord = {}

    def get_positions(self):
        return self._positions

    def get_n_atoms(self):
        self._n_atoms = len(self._positions)
        return self._n_atoms

    def get_symbols(self):
        return self._symbols

    def _add_shape_info(self, shape_label, measure='measure', central_atom=None):
        self._central_atom = central_atom
        self._shape_label = shape_label
        get_measure = 'get_'+measure
        if measure == 'structure':
            self._shape_measures[shape_label]['measure'], self._shape_measures[shape_label][measure] = \
                getattr(shape, get_measure)(self)
        else:
            self._shape_measures[shape_label][measure] = getattr(shape, get_measure)(self)

    def get_shape_measure(self, shape_label, central_atom=None):
        if shape_label not in self._shape_measures:
            self._shape_measures[shape_label] = {}
        if 'measure' not in self._shape_measures[shape_label]:
            self._add_shape_info(shape_label, measure='measure', central_atom=central_atom)
        return self._shape_measures[shape_label]['measure']

    def get_shape_structure(self, shape_label, central_atom=None):
        if shape_label not in self._shape_measures:
            self._shape_measures[shape_label] = {}
        if 'structure' not in self._shape_measures[shape_label]:
            self._add_shape_info(shape_label, measure='structure', central_atom=central_atom)
        return self._shape_measures[shape_label]['structure']

    def get_path_deviation(self, shape_label1, shape_label2, central_atom):
        if shape_label1+'_'+shape_label2 not in self._path_deviation:
            if shape_label2+'_'+shape_label1 not in self._path_deviation:
                labels = shape_label1+'_'+shape_label2
                self._path_deviation[labels] = None
            else:
                labels = shape_label2 + '_' + shape_label1
        else:
            labels = shape_label1 + '_' + shape_label2
        if self._path_deviation[labels] is None:
            Sx = self.get_shape_measure(shape_label1, central_atom)
            Sy = self.get_shape_measure(shape_label2, central_atom)
            self._path_deviation[labels] = shape.get_path_deviation(Sx, Sy, shape_label1, shape_label2)
        return self._path_deviation[labels]

    def get_GenCoord(self, shape_label1, shape_label2, central_atom):
        if shape_label1+'_'+shape_label2 not in self._GenCoord:
            if shape_label2+'_'+shape_label1 not in self._GenCoord:
                labels = shape_label1+'_'+shape_label2
                self._GenCoord[labels] = None
            else:
                labels = shape_label2 + '_' + shape_label1
        else:
            labels = shape_label1 + '_' + shape_label2
        if self._GenCoord[labels] is None:
            Sq = self.get_shape_measure(shape_label1,  central_atom)
            self._GenCoord[labels] = shape.get_GenCoord(Sq, shape_label1, shape_label2)
        return self._GenCoord[labels]
