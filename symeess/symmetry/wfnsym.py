from wfnsympy import WfnSympy
from symeess import tools
import hashlib


class WFNSYM:

    def __init__(self, molecule):

        self._molecule = molecule
        self._Ne_val = self._get_valence_electrons()
        self._wfnsym_dict = {}

        self._results = {}

    def measure(self, label, vector_axis2, vector_axis1=None, center=None):
        if center is None:
            center = [0., 0., 0.]
        if vector_axis1 is None:
            vector_axis1 = [0., 0., 1.]
        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        if hash not in self._results:
            self._basis_to_wfnsym_format()

            self._results[hash] = WfnSympy(NEval=self._Ne_val,
                                           AtLab=self._molecule.geometry.get_symbols(),
                                           shell_type=self._wfnsym_dict['shell_type'],
                                           p_exp=self._wfnsym_dict['p_exponents'],
                                           con_coef=self._wfnsym_dict['con_coefficients'],
                                           p_con_coef=self._wfnsym_dict['p_con_coefficients'],
                                           RAt=self._molecule.geometry.get_positions(),
                                           n_prim=self._wfnsym_dict['n_primitive'],
                                           atom_map=self._wfnsym_dict['atom_map'],
                                           Ca=self._molecule.electronic_structure.coefficients_a,
                                           Cb=self._molecule.electronic_structure.coefficients_b,
                                           RCread=center, VAxis=vector_axis1, VAxis2=vector_axis2,
                                           iCharge=self._molecule.electronic_structure.charge,
                                           iMult=self._molecule.electronic_structure.multiplicity,
                                           group=label.upper(),
                                           do_operation=False)
        return self._results[hash]

    def _get_valence_electrons(self):
        n_valence = 0
        for symbol in self._molecule.geometry.get_symbols():
            n_valence += tools.element_valence_electron(symbol)
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
        for idn, symbol in enumerate(self._molecule.geometry.get_symbols()):
            for orbital_type in self._molecule.electronic_structure.basis[symbol]:
                self._wfnsym_dict['n_primitive'].append(len(orbital_type['p_exponents']))
                self._wfnsym_dict['atom_map'].append(idn+1)
                p_exponents.append(orbital_type['p_exponents'])
                con_coefficients.append(orbital_type['con_coefficients'])
                if len(orbital_type) == 3:
                    p_con_coefficients.append([0. for _ in range(len(orbital_type['con_coefficients']))])
                else:
                    p_con_coefficients.append(orbital_type['p_con_coefficients'])
                for key, kind in typeList.items():
                    if kind[0].upper() == orbital_type['shell_type']:
                        self._wfnsym_dict['shell_type'].append(int(key))
                        break
        self._wfnsym_dict['p_exponents'] = [item for sublist in p_exponents for item in sublist]
        self._wfnsym_dict['con_coefficients'] = [item for sublist in con_coefficients for item in sublist]
        self._wfnsym_dict['p_con_coefficients'] = [item for sublist in p_con_coefficients for item in sublist]