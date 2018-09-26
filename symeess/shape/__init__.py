from symeess.shape import shp, shape_tools
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

    # Function description
    def measure(self, label, central_atom=None):
        c_atom = False
        if central_atom is not None:
            coordinates = shape_tools._order_coordinates(self._coordinates, [central_atom, len(self._coordinates)])
            c_atom = True
        else:
            coordinates = self._coordinates
        if isinstance(label, str):
            reference_structure = shape_tools.get_test_structure(label, central_atom)
        else:
            reference_structure = np.array(label)

        hash = hashlib.md5('{}{}{}'.format(coordinates, c_atom, reference_structure).encode()).hexdigest()
        if hash not in self._measures:
            self._measures[hash] = shp.cshm(coordinates, c_atom, reference_structure)

        return self._measures[hash]

    # Function description
    def structure(self, label, central_atom=None):
        c_atom = False
        if central_atom is not None:
            coordinates = shape_tools._order_coordinates(self._coordinates, [central_atom, len(self._coordinates)])
            c_atom = True
        else:
            coordinates = self._coordinates
        if isinstance(label, str):
            reference_structure = shape_tools.get_test_structure(label, central_atom)
        else:
            reference_structure = np.array(label)

        hash = hashlib.md5('{}{}{}'.format(coordinates, c_atom, reference_structure).encode()).hexdigest()
        if hash not in self._structures:
            self._structures[hash], self._measures[hash] = shp.poly(coordinates, c_atom, reference_structure)

        return self._structures[hash]


