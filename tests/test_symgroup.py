import unittest
from cosymlib import file_io
from numpy import testing
from cosymlib.molecule import Geometry


class TestSymgroupFchk(unittest.TestCase):

    def setUp(self):
        self._structure = file_io.read_generic_structure_file('data/wfnsym/tih4_5d.fchk')
        self._geometry = self._structure.geometry

    def test_symmetry_measure(self):
        print(self._structure.geometry)
        measure = self._geometry.get_symmetry_measure('C3', central_atom=1)
        self.assertAlmostEqual(measure, 0)


class TestSymgroupCycles(unittest.TestCase):

    def setUp(self):
        self._geometry = Geometry(positions=[[ 0.506643354, -1.227657970, 0.000000000],
                                             [ 1.303068499,  0.000000000, 0.000000000],
                                             [ 0.506643354,  1.227657970, 0.000000000],
                                             [-0.926250976,  0.939345948, 0.000000000],
                                             [-0.926250976, -0.939345948, 0.000000000]],
                                  # name='test',
                                  symbols=['C', 'C', 'C', 'C', 'C'],
                                  connectivity_thresh=1.5,
                                  )

    def test_symmetry_measure(self):
        measure = self._geometry.get_symmetry_measure('C5')
        self.assertAlmostEqual(measure, 0.8247502, places=6)
        measure = self._geometry.get_symmetry_measure('C2')
        self.assertAlmostEqual(measure, 0.0, places=6)
        measure = self._geometry.get_symmetry_measure('C3')
        self.assertAlmostEqual(measure, 33.482451, places=6)

    #def test_symmetry_measure_permutation(self):
    #    measure = self._geometry.get_symmetry_measure('C5', fix_permutation=True)
    #    self.assertAlmostEqual(measure, 0.8247502, places=6)

    def test_symmetry_nearest(self):
        nearest = self._geometry.get_symmetry_nearest_structure('C5')
        print(nearest)
        reference = [[ 4.05078542e-01, -1.24670356e+00,  0.00000000e+00],
                     [ 1.31086170e+00, -1.33226763e-16,  0.00000000e+00],
                     [ 4.05078542e-01,  1.24670356e+00,  0.00000000e+00],
                     [-1.06050939e+00,  7.70505174e-01,  0.00000000e+00],
                     [-1.06050939e+00, -7.70505174e-01,  0.00000000e+00]]
        testing.assert_array_almost_equal(nearest, reference, decimal=6)
