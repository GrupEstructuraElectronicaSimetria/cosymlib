from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.simulation import ExtendedHuckel
from warnings import warn


# Gets the parameters defined in the arguments of the function and sets them to Symmetry instance
def set_parameters(func):
    def wrapper(*args, **kwargs):
        args[0]._symmetry.set_parameters(kwargs)
        return func(*args, **kwargs)
    return wrapper


class Molecule:
    def __init__(self, geometry, electronic_structure=None, name=None):

        if not geometry:
            print('No geometry found in the input file, check out input file for possible errors')
            exit()
        if name is None:
            self._name = geometry.name

        self._geometry = geometry
        self._electronic_structure = electronic_structure
        self._symmetry = geometry._symmetry
        self._shape = geometry._shape
        self._symmetry.set_electronic_structure(electronic_structure)

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
            self.symmetry.set_electronic_structure(self._electronic_structure)

        return self._electronic_structure

    # TODO: 'symmetry' and 'shape' properties should be removed. All methods inside these
    # TODO: classes should be mirrored in geometry/molecule class
    @property
    def symmetry(self):
        return self._symmetry

    # TODO: Old method (to be deprecated)
    @set_parameters
    def get_mo_symmetry(self, group, axis=None, axis2=None, center=None):
        return self.symmetry._get_wfnsym_results(group)

    # New ones (to substitute get_mo_symmetry)
    @set_parameters
    def get_mo_irreducible_representations(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.mo_irreducible_representations(group)

    @set_parameters
    def get_wf_irreducible_representations(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.wf_irreducible_representations(group)

    @set_parameters
    def get_mo_overlaps(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.mo_overlaps(group)

    @set_parameters
    def get_wf_overlaps(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.wf_overlaps(group)

    @set_parameters
    def get_symmetry_matrix(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.symmetry_matrix(group)
