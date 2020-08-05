from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.simulation import ExtendedHuckel
from warnings import warn


# Gets the parameters defined in the arguments of the function and sets them to Symmetry instance
def set_parameters(func):
    def wrapper(*args, **kwargs):
        args[0]._symmetry.set_parameters(kwargs)
        args[0]._symmetry.set_electronic_structure(args[0].electronic_structure)
        return func(*args, **kwargs)
    return wrapper


class Molecule:
    def __init__(self, geometry,
                 electronic_structure=None,
                 name=None):

        if not geometry:
            raise Exception('No geometry found in the input file, check out input file for possible errors')
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
            warn('Warning: Electronic structure auto generated from extended Huckel')
            eh = ExtendedHuckel(self.geometry)
            self._electronic_structure = ElectronicStructure(basis=eh.get_basis(),
                                                             orbital_coefficients=[eh.get_mo_coefficients(), []],
                                                             mo_energies=eh.get_mo_energies(),
                                                             multiplicity=eh.get_multiplicity(),
                                                             alpha_electrons=[eh.get_alpha_electrons()],
                                                             beta_electrons=[eh.get_beta_electrons()])
            self.symmetry.set_electronic_structure(self._electronic_structure)

        return self._electronic_structure

    # TODO: 'symmetry' and 'shape' properties should be removed. All methods inside these
    # TODO: classes should be mirrored in geometry/molecule class
    @property
    def symmetry(self):
        return self._symmetry

    # TODO: Old method (to be deprecated)
    @set_parameters
    def OLD_get_mo_symmetry(self, group, axis=None, axis2=None, center=None):
        warn('This method is deprecated', DeprecationWarning)
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

    @set_parameters
    def get_wf_symmetry(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.wf_measure(group)

    @set_parameters
    def get_dens_symmetry(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.dens_measure(group)

    @set_parameters
    def get_symmetry_axes(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.axes(group)

    @set_parameters
    def get_ideal_group_table(self, group, axis=None, axis2=None, center=None):
        return self._symmetry.wf_ideal_group_table(group)
