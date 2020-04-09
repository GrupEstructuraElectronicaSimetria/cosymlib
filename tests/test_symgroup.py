import unittest
from cosymlib import file_io
from numpy import testing


class TestSymgroup(unittest.TestCase):

    def setUp(self):
        self._structure = file_io.read_generic_structure_file('data/wfnsym/tih4_5d.fchk')
        self._geometry = self._structure.geometry

    def test_symmetry_measure(self):
        print(self._structure.geometry)
        measure = self._geometry.get_symmetry_measure('C3', central_atom=1)
        self.assertAlmostEqual(measure, 0)
