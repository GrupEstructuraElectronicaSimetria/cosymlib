from wfnsympy import WfnSympy


class ElectronicStructure:
    def __init__(self, electronic_data, geometry):

        self._wfnsym_dict = {}
        self._geometry = geometry
        key_list = ['Charge', 'Mult', 'N_e', 'shell_type', 'n_primitive', 'atom_map',
                    'p_exponents', 'con_coefficients', 'p_con_coefficients', 'Ca', 'Cb']

        for idn, key in enumerate(key_list[:6]):
            self._wfnsym_dict[key] = [int(j) for j in electronic_data[idn]]

        for idn, key in enumerate(key_list[6:-1]):
            self._wfnsym_dict[key] = [float(j) for j in electronic_data[idn + 6]]

        if not electronic_data[-1]:
            self._wfnsym_dict[key_list[-1]] = [float(j) for j in electronic_data[-2]]

        self._wfnsym_dict['N_Val'] = self._get_valence_electrons()

    def get_wyfsym_measure(self, iGroup, nGroup, VAxis1, VAxis2, RCread):
        if '-2' in self._wfnsym_dict['shell_type']:
            pure_d = True
        else:
            pure_d = False
        results = WfnSympy(Etot=self._wfnsym_dict['N_e'],
                           NEval=self._wfnsym_dict['N_Val'],
                           AtLab=self._geometry.get_symbols(),
                           shell_type=self._wfnsym_dict['shell_type'],
                           p_exp=self._wfnsym_dict['p_exponents'],
                           con_coef=self._wfnsym_dict['con_coefficients'],
                           p_con_coef=self._wfnsym_dict['p_con_coefficients'],
                           RAt=self._geometry.get_positions(),
                           n_prim=self._wfnsym_dict['n_primitive'],
                           atom_map=self._wfnsym_dict['atom_map'],
                           Ca=self._wfnsym_dict['Ca'], Cb=self._wfnsym_dict['Cb'],
                           RCread=RCread, VAxis=VAxis1, VAxis2=VAxis2,
                           iCharge=self._wfnsym_dict['Charge'], iMult=self._wfnsym_dict['Mult'],
                           igroup=iGroup,
                           ngroup=nGroup,
                           do_operation=False,
                           use_pure_d_functions=pure_d)
        # results.print_alpha_mo_IRD()
        # results.print_beta_mo_IRD()
        # results.print_wf_mo_IRD()
        # results.print_CSM()
        # results.print_ideal_group_table()
        # results.print_overlap_mo_alpha()
        # results.print_overlap_mo_beta()
        # results.print_overlap_wf()
        # for i in range(results.dgroup):
        #     results.print_symmetry_operation_matrix(i)
        #     results.print_symmetry_transformed_coordinates(i)
        return results

    def _get_valence_electrons(self):
        n_valence = 0
        for symbol in self._geometry.get_symbols():
            n_valence += atoms_electrons[symbol]
        return n_valence


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