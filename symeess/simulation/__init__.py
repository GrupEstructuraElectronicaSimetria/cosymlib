import symeess
import extended_huckel
import file_io



class ExtendedHuckel:

    def __init__(self, geometry):
        self._EH = extended_huckel.ExtendedHuckel(geometry.get_positions(), geometry.get_symbols(), pure_orbitals=False)

    def get_mo_coefficients(self):
        return self._EH.get_eigenvectors()

    def get_basis(self):
        return self._EH.get_molecular_basis()

    def build_fchk_file(self, name):
        txt_fchk = file_io.build_fchk(self._EH)
        open(name + '.fchk', 'w').write(txt_fchk)


# if __name__ == '__main__':
#     molecules = symeess.file_io.read_input_file('../../examples/coord.xyz')
#     for molecule in molecules:
#         EH = ExtendedHuckel(molecule)
#         txt_fchk = EH.build_fchk_file()
#         open(molecule.get_name() + '.fchk', 'w').write(txt_fchk)
