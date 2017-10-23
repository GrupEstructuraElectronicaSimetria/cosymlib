import unittest
import numpy as np
from symeess import file_io, Symeess


class TestShapeCorFile(unittest.TestCase):

    def test_measure_cor(self):

        molecules = file_io.read('data/coord.cor')
        symeess = Symeess(molecules)

        results = list()
        results.append([molecule.geometry.get_shape_measure('SP-4', central_atom=1) for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('T-4', central_atom=1) for molecule in symeess._molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/results_cor_file')
        self.assertTrue(np.allclose(good_results, calculated_results))

    def test_structure_file(self):
        molecules = file_io.read('data/coord.xyz')
        symeess = Symeess(molecules)

        calculated_results = (symeess._molecules[0].geometry.get_shape_structure('SP-4', central_atom=1))
        good_results = np.loadtxt('data/results_structure')
        self.assertTrue(np.allclose(good_results, calculated_results))

