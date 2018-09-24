from wfnsympy import WfnSympy


class ElectronicStructure:
    def __init__(self,
                 geometry,
                 charge=0,
                 multiplicity=1,
                 basis=None,
                 Ca=None,
                 Cb=None):

        self._wfnsym_dict = {}
        self._charge = charge
        self._mult = multiplicity
        self._basis = basis
        self._Ca = Ca
        if not Cb:
            self._Cb = Ca
        else:
            self._Cb = Cb

        self._geometry = geometry
        # key_list = ['Charge', 'Mult', 'shell_type', 'n_primitive', 'atom_map',
        #             'p_exponents', 'con_coefficients', 'p_con_coefficients', 'Ca', 'Cb']

        self._Ne_val = self._get_valence_electrons()

    def get_wfnsym_measure(self, label, VAxis1, VAxis2, RCread):

        # if -2 in self._wfnsym_dict['shell_type']:
        #     pure_d = True
        # else:
        #     pure_d = False
        self._basis_to_wfnsym_format()

        results = WfnSympy(NEval=self._Ne_val,
                           AtLab=self._geometry.get_symbols(),
                           shell_type=self._wfnsym_dict['shell_type'],
                           p_exp=self._wfnsym_dict['p_exponents'],
                           con_coef=self._wfnsym_dict['con_coefficients'],
                           p_con_coef=self._wfnsym_dict['p_con_coefficients'],
                           RAt=self._geometry.get_positions(),
                           n_prim=self._wfnsym_dict['n_primitive'],
                           atom_map=self._wfnsym_dict['atom_map'],
                           Ca=self._Ca, Cb=self._Cb,
                           RCread=RCread, VAxis=VAxis1, VAxis2=VAxis2,
                           iCharge=self._charge, iMult=self._mult,
                           group=label.upper(),
                           do_operation=False,
                           use_pure_d_functions=False)

        return results

    def _get_valence_electrons(self):
        n_valence = 0
        for symbol in self._geometry.get_symbols():
            n_valence += atoms_electrons[symbol]
        return n_valence

    def _basis_to_wfnsym_format(self):

        typeList = {'0': ['s', 1],
                    '1': ['p', 3],
                    '2': ['d', 6],
                    '3': ['f', 10],
                    '-1': ['sp', 4]}

        self._wfnsym_dict['shell_type'] = []
        self._wfnsym_dict['n_primitive'] = []
        self._wfnsym_dict['atom_map'] = []
        p_exponents = []
        con_coefficients = []
        p_con_coefficients = []
        for idn, symbol in enumerate(self._geometry.get_symbols()):
            for orbital_type in self._basis[str(symbol_map[symbol])]:
                self._wfnsym_dict['n_primitive'].append(len(orbital_type[1]))
                self._wfnsym_dict['atom_map'].append(idn+1)
                p_exponents.append(orbital_type[1])
                con_coefficients.append(orbital_type[2])
                if len(orbital_type) == 3:
                    p_con_coefficients.append([0. for _ in range(len(orbital_type[2]))])
                else:
                    p_con_coefficients.append(orbital_type[3])
                for key, kind in typeList.items():
                    if kind[0].upper() == orbital_type[0]:
                        self._wfnsym_dict['shell_type'].append(int(key))
                        break
        self._wfnsym_dict['p_exponents'] = [item for sublist in p_exponents for item in sublist]
        self._wfnsym_dict['con_coefficients'] = [item for sublist in con_coefficients for item in sublist]
        self._wfnsym_dict['p_con_coefficients'] = [item for sublist in p_con_coefficients for item in sublist]


atoms_electrons = {
        'H': 1,
        'He': 2,
        'Li': 1,
        'Be': 2,
        'B': 3,
        'C': 4,
        'N': 5,
        'O': 6,
        'F': 7,
        'Ne': 8,
        'Na': 1,
        'Mg': 2,
        'Al': 3,
        'Si': 4,
        'P': 5,
        'S': 6,
        'Cl': 7,
        'Ar': 8,
        'K': 1,
        'Ca': 2,
        'Sc': 3,
        'Ti': 4,
        'V': 5,
        'Cr': 6,
        'Mn': 7,
        'Fe': 8,
        'Co': 9,
        'Ni': 10,
        'Cu': 11,
        'Zn': 12,
        'Ga': 13,
        'Ge': 14,
        'As': 15,
        'Se': 16,
        'Br': 17,
        'Kr': 18}

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