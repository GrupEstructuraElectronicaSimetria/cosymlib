from symeess import shape
import numpy as np


class Geometry:
    def __init__(self,
                 symbols=None,
                 positions=None,
                 name=None):

        self._shape_label = 0
        self._central_atom = None
        # self._shape_measures = {}
        self._path_deviation = {}
        self._GenCoord = {}
        self._symbols = []
        self._positions = []

        if name.strip() != '':
            self._name = name
        else:
            self._name = ' ' * 5

        for element in symbols:
            try:
                int(element)
                self._symbols.append(get_element_symbol(int(element)))
            except (ValueError, TypeError):
                self._symbols.append(element.capitalize())
                for ida, a in enumerate(element):
                    try:
                        int(a)
                        self._symbols[-1] = self._symbols[-1][:ida]
                        break
                    except (ValueError, TypeError, IndexError):
                        pass

        try:
            float(positions[0])
            for element in positions:
                self._positions.append(float(element))
            self._positions = list(chunks(self._positions, 3))
        except (ValueError, TypeError, IndexError):
            for element in positions:
                self._positions.append([float(j) for j in element])

        self._positions = np.array(self._positions)
        self._shape = shape.Shape(self._positions)

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

    def get_shape_measure(self, shape_label, central_atom=None):
        return self._shape.measure(shape_label,
                                   central_atom=central_atom)

    def get_shape_structure(self, shape_label, central_atom=None):
        return self._shape.structure(shape_label,
                                     central_atom=central_atom)

    # def get_test_structure(self, shape_label, central_atom=None):
    #     return self._shape.test_structure(shape_label,
    #                                       central_atom=central_atom)

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

    def get_generalized_coordinate(self, shape_label1, shape_label2, central_atom):
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
            self._GenCoord[labels] = shape.get_generalized_coordinate(Sq, shape_label1, shape_label2)
        return self._GenCoord[labels]


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def get_element_symbol(atomic_number):
    for key, Z in symbol_map.items():
        if Z == atomic_number:
            return key


symbol_map = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Uut": 113,
    "Uuq": 114,
    "Uup": 115,
    "Uuh": 116,
    "Uus": 117,
    "Uuo": 118,
}