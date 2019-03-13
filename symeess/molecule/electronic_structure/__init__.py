# from symeess.symmetry import wfnsym
# from symeess import tools
import numpy as np


class ElectronicStructure:
    def __init__(self,
                 charge=0,
                 multiplicity=1,
                 basis=None,
                 orbital_coefficients=None):

        self._charge = charge
        self._multiplicity = multiplicity
        self._basis = basis
        self._Ca = [float(y) for x in orbital_coefficients[0] for y in x]
        if not orbital_coefficients[1]:
            self._Cb = [float(y) for x in orbital_coefficients[0] for y in x]
        else:
            self._Cb = [float(y) for x in orbital_coefficients[1] for y in x]

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
