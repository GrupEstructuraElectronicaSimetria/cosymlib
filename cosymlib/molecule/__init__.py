from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.symmetry.wfnsym import Wfnsym
from cosymlib.symmetry.pointgroup import CalculatePointGroup
from cosymlib.simulation import ExtendedHuckel


class Molecule:

    def __init__(self, geometry, ee=None):

        if not geometry:
            print('No geometry found in the input file, check out input file for possible errors')
            exit()
        self._geometry = geometry
        self._name = geometry.get_name()
        self._electronic_structure = ee
        if ee is not None:
            self._wfnsym = Wfnsym(self)

    def get_name(self):
        return self._name

    # def set_wfnsym(self):
        # self._electronic_structure = ee


    @property
    def geometry(self):
        return self._geometry

    @property
    def electronic_structure(self):
        if self._electronic_structure is None:
            print('No electronic structure found in the input file. Starting a extended-huckel calculation to determine'
                  'the molecular orbital coefficients...')
            eh = ExtendedHuckel(self.geometry)
            self._electronic_structure = ElectronicStructure(basis=eh.get_basis(),
                                                             orbital_coefficients=[eh.get_mo_coefficients(), []],
                                                             mo_energies=eh.get_mo_energies(),
                                                             valence_only=True)
        return self._electronic_structure

    def get_mo_symmetry(self, group, vector_axis1=None, vector_axis2=None, center=None):
        from copy import deepcopy
        wfnsym = deepcopy(Wfnsym(self))
        return wfnsym.results(group, vector_axis1, vector_axis2, center)

    def get_pointgroup(self, tol=0.01):
        return CalculatePointGroup(self._geometry, tolerance=tol).get_point_group()
