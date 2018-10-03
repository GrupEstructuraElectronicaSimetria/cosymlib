from symeess.molecule.geometry import Geometry
from symeess.molecule.electronic_structure import ElectronicStructure
from symmetry.wfnsym import Wfnsym


class Molecule:

    def __init__(self, geometry, ee=None):

        if not geometry:
            print('No geometry found in the input file, check out for possible errors')
            exit()
        self._geometry = geometry
        self._name = geometry.get_name()
        if ee is not None:
            self._electronic_structure = ee
            self._wfnsym = Wfnsym(self)

    def get_name(self):
        return self._name


    @property
    def geometry(self):
        return self._geometry

    @property
    def electronic_structure(self):
        return self._electronic_structure

    def get_mo_symmetry(self, label, vector_axis2, vector_axis1, center):
        return self._wfnsym.get_results(label, vector_axis2, vector_axis1, center)

    # get_shape(*)  ????
    #   return self.geometry.get_shape_measure(*)

    # get_structural_symmetry()  ???


    # def calculate_pointgroup(self):
    #     symmetry.get_pointgroup(self._geometry.get_symbols(), self._geometry.get_positions())