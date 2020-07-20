import huckelpy


class ExtendedHuckel:

    def __init__(self, geometry):
        self._EH = huckelpy.ExtendedHuckel(geometry.get_positions(), geometry.get_symbols())
        self._alpha_electrons = None
        self._beta_electrons = None
        self._total_electrons = self._EH.get_number_of_electrons()

    def get_mo_coefficients(self):
        return self._EH.get_eigenvectors()

    def get_basis(self):
        return self._EH.get_molecular_basis()

    def get_mo_energies(self):
        return self._EH.get_mo_energies()

    def get_multiplicity(self):
        return self._EH.get_multiplicity()

    def get_alpha_electrons(self):
        if self._alpha_electrons is None:
            self._alpha_electrons = self._total_electrons // 2 + self.get_multiplicity() - 1
        return self._alpha_electrons

    def get_beta_electrons(self):
        return self._total_electrons - self._alpha_electrons

    def build_fchk_file(self, name):
        txt_fchk = huckelpy.file_io.build_fchk(self._EH)
        open(name + '.fchk', 'w').write(txt_fchk)
