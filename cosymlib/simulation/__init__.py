import huckelpy


class ExtendedHuckel:

    def __init__(self, geometry):
        self._EH = huckelpy.ExtendedHuckel(geometry.get_positions(), geometry.get_symbols())

    def get_mo_coefficients(self):
        return self._EH.get_eigenvectors()

    def get_basis(self):
        return self._EH.get_molecular_basis()

    def get_mo_energies(self):
        return self._EH.get_mo_energies()

    def build_fchk_file(self, name):
        txt_fchk = huckelpy.file_io.build_fchk(self._EH)
        open(name + '.fchk', 'w').write(txt_fchk)
