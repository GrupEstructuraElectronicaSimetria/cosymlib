import unittest
import numpy as np
from symeess import file_io, Symeess


class TestShapeCorFile(unittest.TestCase):

    def test_example01(self):
        molecules, options = file_io.read('data/shape_examples/example01.xyz')
        symeess = Symeess(molecules)
        results = []
        results.append([molecule.geometry.get_shape_measure('T-4', central_atom=5) for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('SP-4', central_atom=5) for molecule in symeess._molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape_examples/example01_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example02(self):
        molecules, options = file_io.read('data/shape_examples/example02.cor')
        symeess = Symeess(molecules)
        results = []
        results.append([molecule.geometry.get_shape_measure('SP-4', central_atom=1) for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('T-4', central_atom=1) for molecule in symeess._molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape_examples/example02_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example03(self):
        molecules, options = file_io.read('data/shape_examples/example03.xyz', True)
        symeess = Symeess(molecules)
        results = []
        results.append([molecule.geometry.get_shape_measure('T-4', central_atom=1) for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('SP-4', central_atom=1) for molecule in symeess._molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape_examples/example03_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example04(self):
        molecules, options = file_io.read('data/shape_examples/example04.xyz', True)
        symeess = Symeess(molecules)
        results = []
        results.append([molecule.geometry.get_shape_structure('T-4', central_atom=1)
                        for molecule in symeess._molecules])
        calculated_results = np.concatenate((results[0][0], results[0][1]))
        results.append([molecule.geometry.get_shape_structure('SP-4', central_atom=1)
                        for molecule in symeess._molecules])
        calculated_results =np.concatenate((calculated_results , np.concatenate((results[1][0], results[1][1]))))
        results = []
        results.append([molecule.geometry.get_shape_measure('T-4', central_atom=1) for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('SP-4', central_atom=1) for molecule in symeess._molecules])
        calculated_results = [calculated_results, np.column_stack((results[0], results[1]))]
        good_results = np.loadtxt('data/shape_examples/example04_results')
        good_results = [good_results, np.loadtxt('data/shape_examples/example03_results')]
        self.assertTrue(np.allclose(good_results[0], calculated_results[0], atol=1e-3))
        self.assertTrue(np.allclose(good_results[1], calculated_results[1], atol=1e-3))

    def test_example05(self):
        molecules, options = file_io.read('data/shape_examples/example05.xyz', True)
        symeess = Symeess(molecules)
        GenCoord = [molecule.geometry.get_GenCoord('OC-6', 'TPR-6', central_atom=2)
                    for molecule in symeess._molecules]
        results = []
        results.append([molecule.geometry.get_shape_measure('OC-6', central_atom=2)
                        for molecule in symeess._molecules])
        results.append([molecule.geometry.get_shape_measure('TPR-6', central_atom=2)
                        for molecule in symeess._molecules])
        path = [molecule.geometry.get_path_deviation('OC-6', 'TPR-6', central_atom=2)
                    for molecule in symeess._molecules]
        map = symeess.get_shape_map('OC-6', 'TPR-6', central_atom=2, num_points=20)
        good_results = np.loadtxt('data/shape_examples/example05_results_gencoord')
        self.assertTrue(np.allclose(good_results, GenCoord, atol=1e-1))
        good_results = np.loadtxt('data/shape_examples/example05_results_measure')
        calculated_results = np.column_stack((results[0], results[1]))
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-1))
        good_results = np.loadtxt('data/shape_examples/example05_results_devpath')
        self.assertTrue(np.allclose(good_results, path, atol=1e-1))
        good_results = np.loadtxt('data/shape_examples/example05_results_map')
        self.assertTrue(np.allclose(good_results.T, map, atol=1e-1))
