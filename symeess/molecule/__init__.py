from symeess.molecule.geometry import Geometry
from symeess.molecule.electronic_structure import ElectronicStructure
from symeess.symmetry.wfnsym import Wfnsym
from symeess.simulation import ExtendedHuckel


class Molecule:

    def __init__(self, geometry, ee=None):

        if not geometry:
            print('No geometry found in the input file, check out for possible errors')
            exit()
        self._geometry = geometry
        self._name = geometry.get_name()
        self._electronic_structure = ee
        if ee is not None:
            self._wfnsym = Wfnsym(self)

    def get_name(self):
        return self._name

    def set_electronic_structure(self, ee):
        self._electronic_structure = ee
        self._wfnsym = Wfnsym(self)

    @property
    def geometry(self):
        return self._geometry

    @property
    def electronic_structure(self):
        return self._electronic_structure

    def get_mo_symmetry(self, group, vector_axis1, vector_axis2, center):
        if self._electronic_structure is None:
            print('No electronic structure found in the input file. Starting a extended-huckel calculation to determine'
                  'the molecular orbital coefficients...')
            EH = ExtendedHuckel(self.geometry)
            self.set_electronic_structure(ElectronicStructure(basis=EH.get_basis(),
                                                              orbital_coefficients=[EH.get_mo_coefficients(),
                                                                                    []],
                                                              valence_only=True))
        return self._wfnsym.results(group, vector_axis1, vector_axis2, center)

    # get_shape(*)  ????
    #   return self.geometry.get_shape_measure(*)

    # get_structural_symmetry()  ???


    # def calculate_pointgroup(self):
    #     symmetry.get_pointgroup(self._geometry.get_symbols(), self._geometry.get_positions())