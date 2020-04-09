from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.simulation import ExtendedHuckel
from cosymlib.symmetry import Symmetry
from cosymlib.shape import Shape
from warnings import warn


class Molecule:
    def __init__(self, geometry, electronic_structure=None):

        if not geometry:
            print('No geometry found in the input file, check out input file for possible errors')
            exit()
        self._name = geometry.name
        self._geometry = geometry
        self._electronic_structure = electronic_structure
        self._symmetry = None
        self._shape = None

    @property
    def name(self):
        return self._name

    @property
    def geometry(self):
        return self._geometry

    @property
    def electronic_structure(self):
        if self._electronic_structure is None:
            warn('No electronic structure found in the input file.' +
                  'Starting a extended-huckel calculation to determine' +
                  'the molecular orbital coefficients...')
            eh = ExtendedHuckel(self.geometry)
            self._electronic_structure = ElectronicStructure(basis=eh.get_basis(),
                                                             orbital_coefficients=[eh.get_mo_coefficients(), []],
                                                             mo_energies=eh.get_mo_energies(),
                                                             valence_only=True)
        return self._electronic_structure

    @property
    def symmetry(self):
        if self._symmetry is None:
            self._symmetry = Symmetry(self)

        return self._symmetry

    @property
    def shape(self):
        if self._shape is None:
            self._shape = Shape(self)
        return self._shape

    def get_mo_symmetry(self, group, vector_axis1=None, vector_axis2=None, center=None):
        return self.symmetry.get_wfnsym_results(group, vector_axis1, vector_axis2, center)
