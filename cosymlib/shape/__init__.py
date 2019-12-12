from cosymlib.shape import shp, shape_tools
import numpy as np
import hashlib


class Shape:
    def __init__(self, geometry):

        # Allow geometry or molecule to be imported instead of crude Cartesian coordinates
        try:
            self._coordinates = geometry.get_positions()
        except AttributeError:
            try:
                self._coordinates = geometry.geometry.get_positions()
            except AttributeError:
                self._coordinates = geometry

        self._coordinates = np.ascontiguousarray(self._coordinates)

        self._measures = {}
        self._structures = {}
        self._test_structures = {}
        self._path_deviation = {}
        self._gen_coord = {}

    # Function description
    def measure(self, label, central_atom=0, fix_permutation=False):
        hash = hashlib.md5('{}{}'.format(central_atom, label).encode()).hexdigest()
        if hash not in self._measures:
            if isinstance(label, str):
                reference_structure = shape_tools.get_test_structure(label, central_atom)
            else:
                reference_structure = np.array(label)

            if fix_permutation:
                self._measures[hash] = shp.cshm_fix(self._coordinates, reference_structure, central_atom)
            else:
                self._measures[hash] = shp.cshm(self._coordinates, reference_structure, central_atom)

        return self._measures[hash]

    # Function description
    def structure(self, label, central_atom=0, fix_permutation=False):
        hash = hashlib.md5('{}{}'.format(central_atom, label).encode()).hexdigest()
        if hash not in self._structures:
            if isinstance(label, str):
                reference_structure = shape_tools.get_test_structure(label, central_atom)
            else:
                reference_structure = np.array(label)

            if fix_permutation:
                self._structures[hash], self._measures[hash] = shp.poly_fix(self._coordinates, reference_structure,
                                                                            central_atom)
            else:
                self._structures[hash], self._measures[hash] = shp.poly(self._coordinates, reference_structure,
                                                                        central_atom)

        return self._structures[hash]

    def get_path_deviation(self, shape_label1, shape_label2, central_atom=0):

        hash = hashlib.md5('{}{}{}'.format(shape_label1, shape_label2, central_atom).encode()).hexdigest()
        if hash not in self._path_deviation:
            Sx = self.measure(shape_label1, central_atom)
            Sy = self.measure(shape_label2, central_atom)
            new_theta = np.arcsin(np.sqrt(Sx) / 10) + np.arcsin(np.sqrt(Sy) / 10)
            if isinstance(shape_label1, np.ndarray):
                structure_a = shape_label1
            else:
                structure_a = shape_tools.get_test_structure(shape_label1, central_atom=central_atom)
            theta = np.arcsin(np.sqrt(Shape(structure_a).measure(shape_label2, central_atom=len(structure_a))) / 10)
            self._path_deviation[hash] = ((new_theta / theta) - 1) * 100

        return self._path_deviation[hash]

    def get_generalized_coordinate(self, shape_label1, shape_label2, central_atom=0):

        hash = hashlib.md5('{}{}{}'.format(shape_label1, shape_label2, central_atom).encode()).hexdigest()
        if hash not in self._gen_coord:
            Sq = self.measure(shape_label1, central_atom)
            if isinstance(shape_label1, np.ndarray):
                structure_a = shape_label1
            else:
                structure_a = shape_tools.get_test_structure(shape_label1, central_atom=central_atom)
            theta = np.arcsin(np.sqrt(Shape(structure_a).measure(shape_label2, central_atom=len(structure_a))) / 10)
            self._gen_coord[hash] = round(100 * np.arcsin(np.sqrt(Sq) / 10) / theta, 1)

        return self._gen_coord[hash]



