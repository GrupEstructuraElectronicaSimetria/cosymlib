#from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.simulation import ExtendedHuckel
from warnings import warn
from copy import deepcopy


# Gets the parameters defined in the arguments of the function and sets them to Symmetry instance
def set_parameters(func):
    def wrapper(*args, **kwargs):
        args[0]._symmetry.set_parameters(kwargs)
        args[0]._symmetry.set_electronic_structure(args[0].electronic_structure)
        return func(*args, **kwargs)
    return wrapper


class Molecule:
    """
    Main molecule class

    :param geometry: The geometry
    :type geometry: Geometry
    :param electronic_structure: The electronic structure
    :type electronic_structure: ElectronicStructure
    :param name: Molecule name
    :type name: str

    """
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

    def molecule_copy(self):
        return deepcopy(self)

    def set_name(self, name):
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def geometry(self):
        """
        Get the geometry

        :return: The geometry
        :rtype: Geometry
        """
        return self._geometry

    @property
    def electronic_structure(self):
        """
        Get the electronic structure

        :return: The electronic structure
        :rtype: ElectronicStructure
        """
        if self._electronic_structure is None:
            warn('Warning: Electronic structure auto generated from extended Huckel')
            eh = ExtendedHuckel(self.geometry)
            self._electronic_structure = ElectronicStructure(basis=eh.get_basis(),
                                                             orbital_coefficients=[eh.get_mo_coefficients(), []],
                                                             mo_energies=eh.get_mo_energies(),
                                                             multiplicity=eh.get_multiplicity(),
                                                             alpha_occupancy=[1]*eh.get_alpha_electrons(),
                                                             beta_occupancy=[1]*eh.get_beta_electrons())
            self.symmetry.set_electronic_structure(self._electronic_structure)

        return self._electronic_structure

    # TODO: 'symmetry' and 'shape' properties should be removed. All methods inside these
    # TODO: classes should be mirrored in geometry/molecule class
    @property
    def symmetry(self):
        return self._symmetry

    @property
    def shape(self):
        return self._shape

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

    # Mirror methods in geometry
    def get_connectivity(self):
        return self.geometry.get_connectivity()

    def get_positions(self):
        """
        Get the positions in Cartesian coordinates

        :return: the coordinates
        :rtype: list
        """
        return self.geometry.get_positions()

    def get_n_atoms(self):
        """
        Get the number of atoms

        :return: number of atoms
        :rtype: int
        """
        return len(self.get_positions())

    def get_symbols(self):
        """
        Get the atomic elements symbols

        :return: the symbols
        :rtype: list
        """
        return self.geometry.get_symbols()

    def get_pointgroup(self, tol=0.01):
        """
        Get the symmetry point group
        :param tol: The tolerance
        :type tol: float
        :return: The point group
        :rtype: str
        """
        return self.geometry.get_point_group(tol=tol)
