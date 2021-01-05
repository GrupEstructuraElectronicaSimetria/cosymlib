class ElectronicStructure:
    """
    Main Electronic structure class

    :param charge: The charge
    :type charge: int
    :param multiplicity: The multiplicity
    :type multiplicity: int
    :param basis: The basis set
    :type basis: dict
    :param orbital_coefficients: Molecular orbital coefficients
    :type orbital_coefficients: list
    :param mo_energies: Molecular orbital energies
    :type mo_energies: list
    :param alpha_occupancy: Number of alpha electrons
    :type alpha_occupancy: list
    :param beta_occupancy: Number of beta electrons
    :type beta_occupancy: list

    """
    def __init__(self,
                 basis,
                 orbital_coefficients,
                 charge=0,
                 multiplicity=1,
                 mo_energies=None,
                 alpha_occupancy=None,
                 beta_occupancy=None):

        self._charge = charge
        self._multiplicity = multiplicity
        self._basis = basis
        self._Ca = orbital_coefficients[0]
        if len(orbital_coefficients[1]) == 0:
            self._Cb = orbital_coefficients[0]
        else:
            self._Cb = orbital_coefficients[1]

        self._mo_energies = mo_energies

        # self._alpha_occupancy = [1] * int(alpha_occupancy[0])
        # self._beta_occupancy = [1] * int(beta_electrons[0])
        self._alpha_occupancy = alpha_occupancy
        self._beta_occupancy = beta_occupancy
        self._occupancy_consistency()
        # if self._multiplicity > 1:
        #     for i in range(self._multiplicity-1):
        #         self._beta_occupancy.append(0)
        self._total_electrons = sum(self._alpha_occupancy) + sum(self._beta_occupancy)

    def set_occupancy(self, occupancy):
        self.set_alpha_occupancy(occupancy)
        self.set_beta_occupancy(occupancy)

    def set_alpha_occupancy(self, occupancy):
        self._alpha_occupancy = occupancy
        self._occupancy_consistency()
        self._recalculate_charge_multiplicity()

    def set_beta_occupancy(self, occupancy):
        self._beta_occupancy = occupancy
        self._occupancy_consistency()
        self._recalculate_charge_multiplicity()

    def _recalculate_charge_multiplicity(self):
        self._charge = self._total_electrons - (sum(self._alpha_occupancy) + sum(self._beta_occupancy))
        self._update_total_electrons()
        self._multiplicity = abs(sum([a_electron - self._beta_occupancy[ida]
                                      for ida, a_electron in enumerate(self._alpha_occupancy)])) + 1

    def _occupancy_consistency(self):
        if len(self._beta_occupancy) < len(self._alpha_occupancy):
            for i in range(len(self._alpha_occupancy) - len(self._beta_occupancy)):
                self._beta_occupancy.append(0)
        elif len(self._alpha_occupancy) < len(self._beta_occupancy):
            for i in range(len(self._beta_occupancy) - len(self._alpha_occupancy)):
                self._alpha_occupancy.append(0)

    def _update_total_electrons(self):
        self._total_electrons = sum(self._alpha_occupancy) + sum(self._beta_occupancy)

    @property
    def charge(self):
        return self._charge

    @property
    def multiplicity(self):
        return self._multiplicity

    @property
    def basis(self):
        return self._basis

    @property
    def coefficients_a(self):
        return self._Ca

    @property
    def coefficients_b(self):
        return self._Cb

    @property
    def energies(self):
        return self._mo_energies

    @property
    def alpha_occupancy(self):
        return self._alpha_occupancy

    @property
    def beta_occupancy(self):
        return self._beta_occupancy

    @property
    def alpha_electrons(self):
        return sum(self._alpha_occupancy)

    @property
    def beta_electrons(self):
        return sum(self._beta_occupancy)
