#from cosymlib.molecule.geometry import Geometry
from cosymlib.molecule.electronic_structure import ElectronicStructure
from cosymlib.molecule.electronic_structure import ProtoElectronicDensity
from cosymlib.molecule.electronic_structure import ProtoElectronicStructure
from cosymlib.simulation import ExtendedHuckel
from cosymlib.tools import element_to_atomic_number
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
    :type geometry: Geometry, Molecule
    :param electronic_structure: The electronic structure
    :type electronic_structure: ElectronicStructure, str
    :param name: Molecule name
    :type name: str

    """
    def __init__(self, geometry,
                 electronic_structure=None,
                 name=None):

        # handle if geometry is already a molecule
        if isinstance(geometry, Molecule):
            if electronic_structure is None:
                electronic_structure = geometry.electronic_structure
            geometry = geometry.geometry
        if name is None:
            self._name = geometry.name

        self._geometry = geometry
        self._electronic_structure = electronic_structure
        self._symmetry = geometry._symmetry
        self._shape = geometry._shape

    def get_copy(self):
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

        # if None default electronic structure is ExtendedHuckel
        if self._electronic_structure is None:
            self._electronic_structure = 'EH'
            warn('Electronic structure auto generated from Extended-Huckel calculation')

        if self._electronic_structure == 'EH':
            self._electronic_structure = ExtendedHuckel(self.geometry)

        elif self._electronic_structure == 'Dens':
            protodensity = ProtoElectronicDensity(self.geometry)
            self._electronic_structure = ProtoElectronicStructure(basis=protodensity.get_basis(),
                                                                  orbital_coefficients=[protodensity.get_mo_coefficients(),[]])

        return self._electronic_structure

    # TODO: 'symmetry' and 'shape' properties may be removed. All methods inside these
    # TODO: classes may be mirrored in geometry/molecule class? [idea]
    @property
    def symmetry(self):
        return self._symmetry

    @property
    def shape(self):
        return self._shape

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

    def get_charge(self):
        """
        Get the charge of the molecule
        """
        if self._electronic_structure is None:
            return None
        net_electrons = sum([element_to_atomic_number(symbol) for symbol in self.get_symbols()])
        return net_electrons - (sum(self.electronic_structure.alpha_occupancy) +
                                sum(self.electronic_structure.beta_occupancy))
