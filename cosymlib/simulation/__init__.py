from cosymlib.molecule.electronic_structure import ElectronicStructure
from huckelpy.file_io import build_fchk
import huckelpy


class ExtendedHuckel(ElectronicStructure):

    def __init__(self, geometry, charge=0):
        self._EH = huckelpy.ExtendedHuckel(geometry.get_positions(), geometry.get_symbols(), charge=charge)
        alpha_electrons = None
        total_electrons = self._EH.get_number_of_electrons()

        if alpha_electrons is None:
            alpha_electrons = total_electrons // 2 + self._EH.get_multiplicity() - 1

        beta_electrons = total_electrons - alpha_electrons

        super().__init__(basis=self._EH.get_molecular_basis(),
                         orbital_coefficients=[self._EH.get_eigenvectors(), []],
                         alpha_energies=self._EH.get_mo_energies(),
                         beta_energies=[],
                         multiplicity=self._EH.get_multiplicity(),
                         alpha_occupancy=[1] * alpha_electrons,
                         beta_occupancy=[1] * beta_electrons)

    def build_fchk_file(self, name):
        txt_fchk = huckelpy.file_io.build_fchk(self._EH)
        open(name + '.fchk', 'w').write(txt_fchk)
