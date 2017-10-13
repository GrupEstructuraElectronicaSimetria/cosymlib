from molecule.geometry import Geometry


class Molecule:

    def __init__(self, structure=None, ee=None):

        self._geometry = Geometry(structure)

    @property
    def geometry(self):
        return self._geometry

    # @geometry.setter
    # def geometry(self, positions):
    #     self._geometry = Geometry(positions)
    #
    # def get_geometry(self):
    #     return self._geometry.get_positions()


    # molecula.geometry.get_coordinates()
    # melecula.estructure.get_orbital(3)
    #
    # geom = molecule.get_geometry()
    # coordinates = molecule.geometry.get_coorinates()
    # n_atoms = molecule.geometry.get_n_atoms()