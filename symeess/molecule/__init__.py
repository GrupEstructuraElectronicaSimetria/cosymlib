from symeess.molecule.geometry import Geometry


class Molecule:

    def __init__(self, structure_data, electronic_structure=None, name=None):

        self._geometry = Geometry(structure_data)
        self._name = None
        self.set_name(name)

    @property
    def geometry(self):
        return self._geometry

    def set_name(self, name):
        if name.strip() != '':
            self._name = name
        else:
            self._name = ' '*5

    def get_name(self):
        return self._name
