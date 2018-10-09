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
    def measure(self, label, central_atom=0):
        hash = hashlib.md5('{}{}'.format(central_atom, label).encode()).hexdigest()
        if hash not in self._measures:
            if isinstance(label, str):
                reference_structure = shape_tools.get_test_structure(label, central_atom)
            else:
                reference_structure = np.array(label)
                # reference_structure = shape_tools.order_coordinates(reference_structure, [central_atom,
                #                                                                           len(reference_structure)])

            self._measures[hash] = shp.cshm(self._coordinates, reference_structure, central_atom)

        return self._measures[hash]

    # Function description
    def structure(self, label, central_atom=0):
        hash = hashlib.md5('{}{}'.format(central_atom, label).encode()).hexdigest()
        if hash not in self._structures:
            if isinstance(label, str):
                reference_structure = shape_tools.get_test_structure(label, central_atom)
            else:
                reference_structure = np.array(label)
                # reference_structure = shape_tools.order_coordinates(reference_structure, [central_atom,
                #                                                                           len(reference_structure)])

            self._structures[hash], self._measures[hash] = shp.poly(self._coordinates, reference_structure,
                                                                    central_atom)

        return self._structures[hash]


