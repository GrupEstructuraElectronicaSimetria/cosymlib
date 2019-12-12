import unittest
import numpy as np
from cosymlib import Cosymlib
from cosymlib.file_io import old_inputs
from cosymlib import file_io
import cosymlib.shape as shape
import cosymlib.shape.maps as maps


class TestShape(unittest.TestCase):

    def test_example01(self):
        molecules, options = old_inputs.read_old_input('data/shape/example01.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))

        results = []
        results.append([shape.Shape(molecule).measure(reference_polyhedron[0], central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure(reference_polyhedron[1], central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.column_stack((results[0], results[1]))
        nice_measures = [[31.375, 0.97], [33.44, 0.16]]
        self.assertTrue(np.allclose(nice_measures, calculated_results, atol=1e-3))

    def test_example02(self):
        nice_measures = [[5.271, 36.847],
                         [5.184, 36.789],
                         [5.047, 36.698],
                         [5.234, 36.822],
                         [5.1, 36.733]]
        molecules, options = old_inputs.read_old_input('data/shape/example02.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))

        results = [[shape.Shape(molecule).measure(reference_polyhedron[0], options['%central_atom'])
                    for molecule in molecules],
                   [shape.Shape(molecule).measure(reference_polyhedron[1], options['%central_atom'])
                    for molecule in molecules]]
        calculated_results = np.column_stack((results[0], results[1]))
        self.assertTrue(np.allclose(nice_measures, calculated_results, atol=1e-3))

    def test_example03(self):
        nice_measures = [[31.375, 0.97],
                         [33.44, 0.16]]
        molecules, options = old_inputs.read_old_input('data/shape/example03.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))
        results = []
        for reference in reference_polyhedron:
            results.append([shape.Shape(molecule).measure(reference, central_atom=options['%central_atom'])
                            for molecule in molecules])
        calculated_results = np.column_stack((results[0], results[1]))

        self.assertTrue(np.allclose(nice_measures, calculated_results, atol=1e-3))

    def test_example04(self):
        nice_measures = [[31.375, 0.97],
                         [33.44, 0.16]]
        nice_structures = [[2.63702000e+00, 9.00254000e+00, 1.50230800e+01],
                           [4.04611899e+00, 8.76816018e+00, 1.39529277e+01],
                           [2.90730181e+00, 8.13748396e+00, 1.65607229e+01],
                           [1.17328859e+00, 8.36062228e+00, 1.42286293e+01],
                           [2.42137060e+00, 1.07438936e+01, 1.53500401e+01],
                           [-3.47837193e-17, -9.20929468e-18, 2.51183874e-17],
                           [3.77908106e-01, 9.50150647e-01, -8.46343078e-01],
                           [-3.77908106e-01, -9.50150647e-01, 8.46343078e-01],
                           [-1.14966271e+00, 6.33326997e-01, 1.97661211e-01],
                           [1.14966271e+00, -6.33326997e-01, -1.97661211e-01],
                           [2.63702000e+00, 9.00254000e+00, 1.50230800e+01],
                           [3.76728743e+00, 7.46687014e+00, 1.40425715e+01],
                           [4.02414490e+00, 8.97967996e+00, 1.66578721e+01],
                           [1.24989510e+00, 9.02540004e+00, 1.33882879e+01],
                           [1.50675257e+00, 1.05382099e+01, 1.60035885e+01],
                           [9.74091060e-17, 1.01739088e-17, -5.98080786e-17],
                           [5.66862241e-01, 1.42522607e+00, -1.26951447e+00],
                           [-5.66862241e-01, -1.42522607e+00, 1.26951447e+00],
                           [-1.72449417e+00, 9.49990377e-01, 2.96491640e-01],
                           [1.72449417e+00, -9.49990377e-01, -2.96491640e-01]]
        molecules, options = old_inputs.read_old_input('data/shape/example04.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))
        results = []
        results.append([shape.Shape(molecule).structure(reference_polyhedron[0], central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.concatenate((results[0][0], results[0][1]))
        results.append([shape.Shape(molecule).structure(reference_polyhedron[1], central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = np.concatenate((calculated_results, np.concatenate((results[1][0], results[1][1]))))
        results = []
        results.append([shape.Shape(molecule).measure(reference_polyhedron[0], central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure(reference_polyhedron[1], central_atom=options['%central_atom'])
                        for molecule in molecules])
        calculated_results = [calculated_results, np.column_stack((results[0], results[1]))]
        self.assertTrue(np.allclose(nice_structures, calculated_results[0], atol=1e-3))
        self.assertTrue(np.allclose(nice_measures, calculated_results[1], atol=1e-3))

    def test_example05(self):
        nice_gen_coord = [84., 84., 25.4, 79.1, 23.7, 25., 54.6]
        nice_measures = [[12.011, 0.954], [12.012, 0.957], [1.142, 11.245], [10.707, 1.236], [0.993, 12.826],
                         [1.105, 11.783], [5.203, 5.085]]
        nice_dev_path = [7.2, 7.2, 6.5, 5.5, 10.6, 8.1, 8.6]
        nice_path_coordinates = np.array([[0., 16.737],
                                          [0.035, 15.355],
                                          [0.144, 14.002],
                                          [0.329, 12.681],
                                          [0.593, 11.398],
                                          [0.94, 10.158],
                                          [1.371, 8.966],
                                          [1.89, 7.828],
                                          [2.497, 6.748],
                                          [3.195, 5.732],
                                          [3.984, 4.785],
                                          [4.865, 3.911],
                                          [5.837, 3.116],
                                          [6.901, 2.403],
                                          [8.055, 1.777],
                                          [9.298, 1.241],
                                          [10.626, 0.798],
                                          [12.038, 0.45],
                                          [13.53, 0.201],
                                          [15.097, 0.05],
                                          [16.737, 0.]])

        molecules, options = old_inputs.read_old_input('data/shape/example05.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))

        results = []
        results.append([shape.Shape(molecule).measure(reference_polyhedron[0], central_atom=options['%central_atom'])
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure(reference_polyhedron[1], central_atom=options['%central_atom'])
                        for molecule in molecules])
        dev_path = [molecule.get_path_deviation(reference_polyhedron[0], reference_polyhedron[1],
                                                central_atom=options['%central_atom'])
                    for molecule in molecules]
        gen_coord = [molecule.get_generalized_coordinate(reference_polyhedron[0], reference_polyhedron[1],
                                                         central_atom=options['%central_atom'])
                     for molecule in molecules]
        map = maps.get_shape_map(reference_polyhedron[0], reference_polyhedron[1], num_points=20)

        self.assertTrue(np.allclose(nice_gen_coord, gen_coord, atol=1e-1))
        calculated_results = np.column_stack((results[0], results[1]))
        self.assertTrue(np.allclose(nice_measures, calculated_results, atol=1e-1))
        self.assertTrue(np.allclose(nice_dev_path, dev_path, atol=1e-1))
        self.assertTrue(np.allclose(nice_path_coordinates.T, map[:2], atol=1e-1))

    def test_example06(self):
        nice_dev_path = [0.9, 3.7, 4.1, 3.7, 3.1, 0.1, 0.7, 3.9, 4., 1.7, 1.9, 2.6, 1.4,
                         2.1, 1.5, 1.4, 4.3, 3.4, 1.1, 0.9, 1.6, 1.6, 0.1, 1.5, 2.5, 0.2,
                         0.6, 3.1, 1.6, 1.1, 0.5, 1.1, 0.5, 1.6, 0.3, 4.8, 3.4, 3., 4.4,
                         0.6, 0.8, 1., 0.6, 2.8, 3.3, 0.5, 0.6, 2.1, 3.5, 2.6, 4.7, 0.6,
                         0.7, 0.8, 4., 4.8, 1.3, 1.5, 1.3, 0.9, 2.2, 2.5, 4.6, 4.7, 4.7,
                         2.4, 0.8, 0., 1.3, 1.4, 1.4, 0.5, 1.7, 3.4, 4., 2.2, 0.4, 0.,
                         1.5, 1.8, 1., 4.5, 1.5, 1.5, 1.6, 0.7, 0.6, 2.3, 0.3, 0.2, 0.3,
                         4.2, 0.2, 0.4, 0.9, 0.3, 0.4, 0., 3.2, 1.4, 1.8, 1., 1.3, 0.6,
                         2., 0.3, 0., 0.1, 0.6, 2.1, 4.2, 1.6, 2.8, 1.5, 2.6, 2.8, 1.7,
                         0.3, 0.7, 0.5, 0.6, 0.6, 2.5, 1.5, 0.8, 0.3, 0.2, 0.8, 0.7, 3.6,
                         3.8, 4.5, 2.8, 1.5, 3.8, 1.8, 3.9, 1.5, 0.2, 0.4, 0.8, 2.5, 0.5,
                         2.6, 1.3, 0.8, 1.2, 1.9, 3.9, 2.4, 0.9, 0.4, 2.6, 0., 2.4, 1.7,
                         4.4, 1.6, 3.5, 3.8, 3.6, 4.7, 4.8, 3.7, 1.8, 1.6, 2.3, 3.7, 3.,
                         1.3, 0.9, 2.8, 0.9, 4.1, 3.2, 4.6, 1.7, 0.8, 1., 0.7, 0.6, 0.2,
                         4.6, 1.5, 1.5, 1.7, 2.5, 2.9, 0.5, 0.8, 1.3]
        molecules, options = old_inputs.read_old_input('data/shape/example06.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))

        symobj = Cosymlib(molecules)
        shape_measure, devpath, GenCoord = symobj.get_path_parameters(reference_polyhedron[0], reference_polyhedron[1],
                                                                      central_atom=options['%central_atom'], maxdev=5.0)
        self.assertTrue(np.allclose(nice_dev_path, devpath, atol=1e-1))

    def test_example07(self):
        molecules, options = old_inputs.read_old_input('data/shape/example06.dat')
        reference_polyhedron = []
        for number in options['%labels']:
            reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))
        symobj = Cosymlib(molecules)
        shape_measure, devpath, GenCoord = symobj.get_path_parameters(reference_polyhedron[0], reference_polyhedron[1],
                                                                      central_atom=options['%central_atom'], maxdev=10.0,
                                                                      maxgco=60, mingco=40)
        self.assertTrue(np.allclose(0.2, devpath, atol=1e-1))

    def test_example08(self):
        nice_measures = [[1.976, 3.699],
                         [6.955, 0.602]]
        molecules, options = old_inputs.read_old_input('data/shape/example08.dat')
        central_atom = options['%central_atom']
        reference_polyhedron = []
        if options['%labels'] != 0:
            for number in options['%labels']:
                if int(number) == 0:
                    for ref in file_io.get_molecule_from_file_ref('data/shape/example08.ref'):
                        reference_polyhedron.append(ref)
                else:
                    reference_polyhedron.append(shape.shape_tools.get_shape_label(int(number), options['%n_atoms']))

        results = []
        results.append([shape.Shape(molecule).measure(reference_polyhedron[0], central_atom=central_atom)
                        for molecule in molecules])
        results.append([shape.Shape(molecule).measure(reference_polyhedron[1].get_positions(),
                                                      central_atom=central_atom) for molecule in molecules])

        calculated_results = np.column_stack((results[0], results[1]))
        self.assertTrue(np.allclose(nice_measures, calculated_results, atol=1e-2))

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
