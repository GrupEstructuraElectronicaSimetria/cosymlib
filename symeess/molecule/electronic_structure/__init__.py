class ElectronicStructure:
    def __init__(self,
                 charge=0,
                 multiplicity=1,
                 basis=None,
                 orbital_coefficients=None,
                 mo_energies=None,
                 valence_only=False):

        self._charge = charge
        self._multiplicity = multiplicity
        self._basis = basis
        self._valence_only = valence_only
        self._Ca = orbital_coefficients[0]
        if not orbital_coefficients[1]:
            self._Cb = None
        else:
            self._Cb = orbital_coefficients[1]
        self._mo_energies = mo_energies

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
    def valence_only(self):
        return self._valence_only
