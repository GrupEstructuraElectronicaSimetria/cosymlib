from cosymlib.shape import shp, tools
import numpy as np


def _get_key(central_atom, reference, fix_permutation=False, reference_2=''):
    if isinstance(reference, str):
        label_key = reference.lower()
    else:
        label_key = np.array2string(reference.get_positions(), precision=10)

    if isinstance(reference_2, str):
        label2_key = reference_2.lower()
    else:
        label2_key = np.array2string(reference_2.get_positions(), precision=10)

    central_atom_key = int(central_atom)
    fix_permutation_key = str(bool(fix_permutation))

    return central_atom_key, label_key, label2_key, fix_permutation_key


class Shape:
    """
    Shape main class

    :param structure: a geometry, molecule or array type object
    :type structure: Geometry, Molecule, np.array
    """
    def __init__(self, structure):

        # Allow geometry or molecule to be imported instead of crude Cartesian coordinates
        try:
            self._coordinates = structure.get_positions()
        except AttributeError:
            # Try to get from numpy array
            self._coordinates = structure

        self._coordinates = np.ascontiguousarray(self._coordinates)

        self._measures = {}
        self._structures = {}
        self._test_structures = {}
        self._path_deviation = {}
        self._gen_coord = {}

    def measure(self, reference, central_atom=0, fix_permutation=False):
        """
        Get shape measure

        :param reference: Reference shape label or Geometry object
        :type reference: str
        :param central_atom: Central atom position
        :type central_atom: int
        :param fix_permutation: Do not permute atoms during shape calculations
        :type fix_permutation: bool
        :return: The measure
        :rtype: float
        """
        key = _get_key(central_atom, reference, fix_permutation=fix_permutation)
        if key not in self._measures:
            if isinstance(reference, str):
                reference_structure = tools.get_reference_structure(reference, central_atom)
            else:
                reference_structure = reference
            reference_coordinates = reference_structure.get_positions()

            if len(self._coordinates) != len(reference_coordinates):
                raise Exception('Reference and input structures have different number of atoms')

            if fix_permutation:
                self._measures[key] = shp.cshm_fix(self._coordinates, reference_coordinates, central_atom)
            else:
                self._measures[key] = shp.cshm(self._coordinates, reference_coordinates, central_atom)

        return self._measures[key]

    def structure(self, reference, central_atom=0, fix_permutation=False):
        """
        Get the nearest structure to reference

        :param reference: Reference shape label or Geometry object
        :type reference: str, Geometry
        :param central_atom: Central atom position
        :type central_atom: int, Geometry
        :param fix_permutation: Do not permute atoms during shape calculations
        :type fix_permutation: bool
        :return: The structure
        :rtype: Structure
        """
        key = _get_key(central_atom, reference, fix_permutation=fix_permutation)
        if key not in self._structures:
            if isinstance(reference, str):
                reference_structure = tools.get_reference_structure(reference, central_atom)
            else:
                reference_structure = reference

            reference_coordinates = reference_structure.get_positions()

            if len(self._coordinates) != len(reference_coordinates):
                raise Exception('Reference and input structures have different number of atoms')

            if fix_permutation:
                self._structures[key], self._measures[key] = shp.poly_fix(self._coordinates,
                                                                          reference_coordinates,
                                                                          central_atom)
            else:
                self._structures[key], self._measures[key] = shp.poly(self._coordinates,
                                                                      reference_coordinates,
                                                                      central_atom)

        return self._structures[key]

    def path_deviation(self, shape_label1, shape_label2, central_atom=0):
        """
        Get the path deviation

        :param shape_label1: First shape reference label or Geometry object
        :type shape_label1: str, Geometry
        :param shape_label2: Second shape reference label or Geometry object
        :type shape_label2: str, Geometry
        :param central_atom: Central atom position
        :type central_atom: int
        :return: The path deviation
        :rtype: float
        """
        # TODO: Someone improve the description of this function

        key = _get_key(central_atom, shape_label1, reference_2=shape_label2)
        if key not in self._path_deviation:
            Sx = self.measure(shape_label1, central_atom)
            Sy = self.measure(shape_label2, central_atom)
            new_theta = np.arcsin(np.sqrt(Sx) / 10) + np.arcsin(np.sqrt(Sy) / 10)
            if isinstance(shape_label1, str):
                structure_a = tools.get_reference_structure(shape_label1, central_atom=central_atom)
            else:
                structure_a = shape_label1

            theta = np.arcsin(np.sqrt(Shape(structure_a).measure(shape_label2, central_atom=structure_a.get_n_atoms())) / 10)
            self._path_deviation[key] = ((new_theta / theta) - 1) * 100

        return self._path_deviation[key]

    def generalized_coordinate(self, shape_label1, shape_label2, central_atom=0):
        """
        Get the generalized coordinate

        :param shape_label1: First shape reference label or Geometry object
        :type shape_label1: str, Geometry
        :param shape_label2: Second shape reference label or Geometry object
        :type shape_label2: str, Geometry
        :param central_atom: Central atom position
        :type central_atom: int
        :return: The generalized coordinate
        """

        key = _get_key(central_atom, shape_label1, reference_2=shape_label2)
        if key not in self._gen_coord:
            Sq = self.measure(shape_label1, central_atom)
            if isinstance(shape_label1, str):
                structure_a = tools.get_reference_structure(shape_label1, central_atom=central_atom)
            else:
                structure_a = shape_label1

            theta = np.arcsin(np.sqrt(Shape(structure_a).measure(shape_label2, central_atom=structure_a.get_n_atoms())) / 10)
            self._gen_coord[key] = round(100 * np.arcsin(np.sqrt(Sq) / 10) / theta, 1)

        return self._gen_coord[key]
