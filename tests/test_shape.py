import unittest
import numpy as np
from cosym import Cosymlib
from cosym.file_io import old_inputs
import cosym.shape as shape
import cosym.shape.maps as maps


class TestShapeCorFile(unittest.TestCase):

    def test_example01(self):
        molecules, options = old_inputs.read_old_input('data/shape/example01.dat')
        results = []
        results.append([shape.Shape(molecule).measure('T-4', central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure('SP-4', central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape/example01_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example02(self):
        molecules, options = old_inputs.read_old_input('data/shape/example02.dat')
        results = []
        results.append([shape.Shape(molecule).measure('SP-4', options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure('T-4', options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape/example02_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example03(self):
        molecules, options = old_inputs.read_old_input('data/shape/example03.dat')
        results = []
        for number in options['%labels']:
            label = shape.shape_tools.get_shape_label(int(number), options['%n_atoms'])
            results.append([shape.Shape(molecule).measure(label, central_atom=options['%central_atom'])
                            for molecule in molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        good_results = np.loadtxt('data/shape/example03_results')
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-3))

    def test_example04(self):
        molecules, options = old_inputs.read_old_input('data/shape/example04.dat')
        results = []
        results.append([shape.Shape(molecule).structure('T-4',central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.concatenate((results[0][0], results[0][1]))
        results.append([shape.Shape(molecule).structure('SP-4', central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.concatenate((calculated_results , np.concatenate((results[1][0], results[1][1]))))
        results = []
        results.append([shape.Shape(molecule).measure('T-4', central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure('SP-4', central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = [calculated_results, np.column_stack((results[0], results[1]))]
        good_results = np.loadtxt('data/shape/example04_results')
        good_results = [good_results, np.loadtxt('data/shape/example03_results')]
        self.assertTrue(np.allclose(good_results[0], calculated_results[0], atol=1e-3))
        self.assertTrue(np.allclose(good_results[1], calculated_results[1], atol=1e-3))

    def test_example05(self):
        molecules, options = old_inputs.read_old_input('data/shape/example05.dat')


        results = []
        results.append([shape.Shape(molecule).measure('OC-6', central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure('TPR-6', central_atom=options['%central_atom'])
                        for molecule in molecules])
        path = [molecule.get_path_deviation('OC-6', 'TPR-6', central_atom=options['%central_atom'])
                for molecule in molecules]
        GenCoord = [molecule.get_generalized_coordinate('OC-6', 'TPR-6', central_atom=options['%central_atom'])
                    for molecule in molecules]
        map = maps.get_shape_map('OC-6', 'TPR-6', num_points=20)

        good_results = np.loadtxt('data/shape/example05_results_gencoord')
        self.assertTrue(np.allclose(good_results, GenCoord, atol=1e-1))
        good_results = np.loadtxt('data/shape/example05_results_measure')
        calculated_results = np.column_stack((results[0], results[1]))
        self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-1))
        good_results = np.loadtxt('data/shape/example05_results_devpath')
        self.assertTrue(np.allclose(good_results, path, atol=1e-1))
        good_results = np.loadtxt('data/shape/example05_results_map')
        self.assertTrue(np.allclose(good_results.T, map, atol=1e-1))

# Cal esperar a tenir els inputs sencers dels exemples ja que ara no els tinc a ma
    def test_example06(self):
        molecules, options = old_inputs.read_old_input('data/shape/example06.dat')

        symobj = Cosymlib(molecules)
        shape, devpath, GenCoord = symobj.get_path_parameters('SP-4', 'T-4',
                                                               central_atom=options['%central_atom'], maxdev=5.0)
        good_results = np.loadtxt('data/shape/example06_results')
        self.assertTrue(np.allclose(good_results, devpath, atol=1e-1))

    def test_example07(self):
        molecules, options = old_inputs.read_old_input('data/shape/example06.dat')
        symobj = Cosymlib(molecules)
        shape, devpath, GenCoord = symobj.get_path_parameters('SP-4', 'T-4',
                                                               central_atom=options['%central_atom'], maxdev=10.0,
                                                               maxgco=60, mingco=40)
        good_results = np.loadtxt('data/shape/example07_results')
        self.assertTrue(np.allclose(good_results, devpath, atol=1e-1))

    # def test_example08(self):
    #     molecules, options = file_io.read_old_input('data/shape/example08.dat')
    #     custom_ref_structure = file_io.read_input_file('data/shape/example08.ref')
    #
    #     results = []
    #     results.append([shape.Shape(molecule).measure('TPR-6', central_atom=int(options[0][1]))
    #                     for molecule in molecules])
    #     results.append([shape.Shape(molecule).measure(custom_ref_structure[0].get_positions(),
    #                                                   central_atom=int(options[0][1])) for molecule in molecules])
    #
    #     calculated_results = np.column_stack((results[0], results[1]))
    #     good_results = np.loadtxt('data/shape/example08_results_measure')
    #     print(calculated_results)
    #     print(good_results)
    #     self.assertTrue(np.allclose(good_results, calculated_results, atol=1e-2))

    # De l'exemple 9-12 tots necessiten la comanda fixperm que no esta inclosa encara en el programa
    # def test_example09(self):
    #     # fixperm
    #     pass

    # def test_example13(self):
    #
    #     central_atom = 1
    #     ref_str1 = shape.shape_tools.get_test_structure('EP-9', central_atom)
    #     ref_str1 = shape.shape_tools.order_coordinates(ref_str1, [len(ref_str1), central_atom])
    #     ref_str2 = shape.shape_tools.get_test_structure('CSAPR-9', central_atom)
    #     ref_str2 = shape.shape_tools.order_coordinates(ref_str2, [len(ref_str2), central_atom])
    #     good_results = np.loadtxt('data/shape/example09_results')
    #     self.assertTrue(np.allclose(good_results, ref_str1, atol=1e-3))

