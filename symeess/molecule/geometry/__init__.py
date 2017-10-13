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
        self._shape_references = {}

    def get_positions(self):
        return self._positions

    def get_n_atoms(self):
        self._n_atoms = len(self._positions)
        return self._n_atoms

    def get_symbols(self):
        return self._symbols

    def set_measure(self, shape_ideal, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_ideal
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_references[shape_ideal]['measure'] = measure(self.get_positions(), n_atoms,
                                                                 self._shape_ideal, self._central_atom)

    def set_ideal_structure(self, shape_ideal, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_ideal
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_references[shape_ideal]['measure'], self._shape_references[shape_ideal]['structure'] = \
            structure_measure(self.get_positions(), n_atoms,self._shape_ideal, self._central_atom)

    def set_test_structure(self, shape_ideal, central_atom):
        self._central_atom = central_atom
        self._shape_ideal = shape_ideal
        n_atoms = self.get_n_atoms()
        if self._central_atom:
            n_atoms = self.get_n_atoms() - 1
        self._shape_test_structure = test_structure(self.get_positions(), n_atoms,
                                                    self._shape_ideal, self._central_atom)

    def get_measure(self, shape_reference, central_atom):
        if shape_reference not in self._shape_references:
            self._shape_references[shape_reference] = {}
        if 'measure' not in self._shape_references[shape_reference]:
            self.set_measure(shape_reference, central_atom)
        return self._shape_references[shape_reference]['measure']

    def get_ideal_structure(self, shape_reference, central_atom):
        if shape_reference not in self._shape_references:
            self._shape_references[shape_reference] = {}
        if 'structure' not in self._shape_references[shape_reference]:
            self.set_ideal_structure(shape_reference, central_atom)
        return self._shape_references[shape_reference]['structure']

    def get_test_structure(self, shape_reference, central_atom):
        self.set_test_structure(shape_reference, central_atom)
        return self._shape_test_structure


# shape_references = {'L-2': [None, None], 'vT-2': [None, None], 'vOC-2': [None, None],
#                     'TP-3': [None, None], 'vT-3': [None, None], 'fvOC-3': [None, None],
#                     'mvOC-3': [None, None], 'SP-4': [None, None], 'T-4': [None, None],
#                     'SS-4': [None, None], 'PP-5': [None, None], 'vOC-5': [None, None],
#                     'TBPY-5': [None, None], 'SPY-5': [None, None], 'JTBPY-5': [None, None],
#                     'HP-6': [None, None], 'PPY-6': [None, None],
#                     'OC-6': [None, None], 'TPR-6': [None, None], 'JPPY-6': [None, None],
#                     'HP-7': [None, None], 'HPY-7': [None, None], 'PBPY-7': [None, None],
#                     'COC-7': [None, None], 'CTPR-7': [None, None], 'JPBPY-7': [None, None],
#                     'JETPY-7': [None, None], 'OP-8': [None, None], 'HPY-8': [None, None],
#                     'HBPY-8': [None, None], 'CU-8': [None, None], 'SAPR-8': [None, None],
#                     'TDD-8': [None, None], 'JGBF-8': [None, None], 'JETBPY-8': [None, None],
#                     'JBTPR-8': [None, None], 'BTPR-8': [None, None], 'JSD-8': [None, None],
#                     'TT-8': [None, None], 'ETBPY-8': [None, None], 'EP-9': [None, None],
#                     'OPY-9': [None, None], 'HBPY-9': [None, None], 'JTC-9': [None, None],
#                     'JCCU-9': [None, None], 'CCU-9': [None, None], 'JCSAPR-9': [None, None],
#                     'CSAPR-9': [None, None], 'JTCTPR-9': [None, None], 'TCTPR-9': [None, None],
#                     'JTDIC-9': [None, None], 'HH-9': [None, None], 'MFF-9': [None, None],
#                     'DP-10': [None, None], 'EPY-10': [None, None], 'OBPY-10': [None, None],
#                     'PPR-10': [None, None], 'PAPR-10': [None, None], 'JBCCU-10': [None, None],
#                     'JBCSAPR-10': [None, None], 'JMBIC-10': [None, None], 'JATDI-10': [None, None],
#                     'JSPC-10': [None, None], 'SDD-10': [None, None], 'TD-10': [None, None],
#                     'HD-10': [None, None], 'HP-11': [None, None], 'DPY-11': [None, None],
#                     'EBPY-11': [None, None], 'JCPPR-11': [None, None], 'JCPAPR-11': [None, None],
#                     'JAPPR-11': [None, None], 'JASPC-11': [None, None], 'DP-12': [None, None],
#                     'HPY-12': [None, None], 'DBPY-12': [None, None], 'HPR-12': [None, None],
#                     'HAPR-12': [None, None], 'TT-12': [None, None], 'COC-12': [None, None],
#                     'ACOC-12': [None, None], 'IC-12': [None, None], 'JSC-12': [None, None],
#                     'JEPBPY-12': [None, None], 'JBAPPR-12': [None, None], 'JSPMC-12': [None, None],
#                     'DD-20': [None, None], 'TCU-24': [None, None], 'TOC-24': [None, None]}
