from molecule.geometry import Geometry


class Molecule:

    def __init__(self, structure=None, ee=None, name=None):

        self._geometry = Geometry(structure)
        self._name = ''
        self.set_name(name)

    @property
    def geometry(self):
        return self._geometry

    def set_name(self, name):
        if name is not None:
            self._name = name

    def get_name(self):
        return self._name