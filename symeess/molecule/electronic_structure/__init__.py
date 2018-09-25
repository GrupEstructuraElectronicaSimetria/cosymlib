from wfnsympy import WfnSympy
from symeess import tools


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
        self._Ne_val = self._get_valence_electrons()

    def get_wfnsym_measure(self, label, VAxis1, VAxis2, RCread):

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
        for idn, symbol in enumerate(self._geometry.get_symbols()):
            for orbital_type in self._basis[symbol]:
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
