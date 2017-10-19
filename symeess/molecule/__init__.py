from molecule.geometry import Geometry


class Molecule:

    def __init__(self, structure=None, ee=None):

        self._geometry = Geometry(structure)

    @property
    def geometry(self):
        return self._geometry
