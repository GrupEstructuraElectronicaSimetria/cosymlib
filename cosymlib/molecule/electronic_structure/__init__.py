class ElectronicStructure:
    def __init__(self,
                 charge=0,
                 multiplicity=1,
                 basis=None,
                 orbital_coefficients=None,
                 mo_energies=None,
                 alpha_electrons=None,
                 beta_electrons=None):

        self._charge = charge
        self._multiplicity = multiplicity
        self._basis = basis
        self._Ca = orbital_coefficients[0]
        if len(orbital_coefficients[1]) == 0:
            self._Cb = orbital_coefficients[0]
        else:
            self._Cb = orbital_coefficients[1]

        self._mo_energies = mo_energies

        self._alpha_occupancy = [1] * int(alpha_electrons[0])
        self._beta_occupancy = [1] * int(beta_electrons[0])
        if self._multiplicity > 1:
            for i in range(self._multiplicity-1):
                self._beta_occupancy.append(0)

    def set_alpha_occupancy(self, occupancy, restricted=False):
        self._alpha_occupancy = occupancy
        if restricted:
            self.set_beta_occupancy(occupancy)

    def set_beta_occupancy(self, occupancy):
        self._beta_occupancy = occupancy

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

    # @property
    # def valence_only(self):
    #     return self._valence_only

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
