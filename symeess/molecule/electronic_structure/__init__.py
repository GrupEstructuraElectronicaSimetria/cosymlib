# from symeess.symmetry import wfnsym
# from symeess import tools


class ElectronicStructure:
    def __init__(self,
                 charge=0,
                 multiplicity=1,
                 basis=None,
                 orbital_coefficients=None):

        self._charge = charge
        self._mult = multiplicity
        self._basis = basis
        self._Ca = [float(i) for i in orbital_coefficients[0]]
        if not orbital_coefficients[1]:
            self._Cb = [float(i) for i in orbital_coefficients[0]]
        else:
            self._Cb = [float(i) for i in orbital_coefficients[1]]

    @property
    def charge(self):
        return self._charge

    @property
    def multiplicity(self):
        return self._mult

    @property
    def basis(self):
        return self._basis

    @property
    def coefficients_a(self):
        return self._Ca

    @property
    def coefficients_b(self):
        return self._Cb
