from symeess.molecule.geometry import Geometry
from symeess.molecule.electronic_structure import ElectronicStructure
# import symmetry


class Molecule:

    def __init__(self, structure_data, electronic_structure=None, name=None):
        if not structure_data:
            print('No molecule found in the input file, check out for possible errors')
            quit()
        self._geometry = Geometry(structure_data)
        self._name = None
        self.set_name(name)
        if electronic_structure is not None:
            self._electronic_structure = ElectronicStructure(electronic_structure, self.geometry)

    def set_name(self, name):
        if name.strip() != '':
            self._name = name
        else:
            self._name = ' '*5

    def get_name(self):
        return self._name


    @property
    def geometry(self):
        return self._geometry

    @property
    def electronic_structure(self):
        return self._electronic_structure

    # def calculate_pointgroup(self):
    #     symmetry.get_pointgroup(self._geometry.get_symbols(), self._geometry.get_positions())