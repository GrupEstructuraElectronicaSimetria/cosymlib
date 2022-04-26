import unittest
from cosymlib import file_io
from numpy import testing
import os


data_dir = os.path.join(os.path.dirname(__file__), 'data')


class TestPointGroup(unittest.TestCase):

    def test_reading_xyz(self):
        methane = [[ 0.0000,  0.0000,  0.0000],
                   [ 0.5288,  0.1610,  0.9359],
                   [ 0.2051,  0.8240, -0.6786],
                   [ 0.3345, -0.9314, -0.4496],
                   [-1.0685, -0.0537,  0.1921]]
        ammonium = [[ 0.0000,  0.0000,  0.1000],
                    [ 0.5288,  0.1610,  0.9359],
                    [ 0.2051,  0.8240, -0.6786],
                    [ 0.3345, -0.9314, -0.4496],
                    [-1.0685, -0.0537,  0.1921]]
        structure = file_io.get_geometry_from_file_xyz(data_dir + '/file_io/test.xyz', read_multiple=True)
        testing.assert_array_equal(structure[0].get_positions(), methane)
        testing.assert_array_equal(structure[1].get_positions(), ammonium)
        self.assertEqual(structure[0].get_symbols(), ['C', 'H', 'H', 'H', 'H'])
        self.assertEqual(structure[1].get_symbols(), ['N', 'H', 'H', 'H', 'H'])
        self.assertEqual(structure[0].get_n_atoms(), 5)
        self.assertEqual(structure[1].get_n_atoms(), 5)

    def test_reading_cor(self):
        structure1 = [[16.04450,  0.00000,  0.00000],
                      [17.11627, -1.46619, -0.92483],
                      [17.11627,  1.46619,  0.92483],
                      [14.97273, -1.46619, -0.92483],
                      [14.97273,  1.46619,  0.92483]]
        structure2 = [[6.35325,  1.50775,  0.00000],
                      [4.67365,  2.63675,  0.26075],
                      [7.85296,  1.80990,  1.32660],
                      [4.85354,  1.20560, -1.32660],
                      [8.03285,  0.37875, -0.26075]]
        structures = file_io.get_geometry_from_file_cor(data_dir +'/file_io/test.cor', read_multiple=True)
        testing.assert_array_equal(structures[0].get_positions(), structure1)
        testing.assert_array_equal(structures[1].get_positions(), structure2)
        self.assertEqual(structures[0].get_symbols(), ['Pd', 'N', 'N', 'N', 'N'])
        self.assertEqual(structures[1].get_symbols(), ['Pt', 'N', 'N', 'N', 'N'])
        self.assertEqual(structures[0].get_n_atoms(), 5)
        self.assertEqual(structures[1].get_n_atoms(), 5)

    def test_reading_fchk(self):
        alpha_occupancy = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        beta_occupancy = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        alpha_electrons = 18
        beta_electrons = 18
        molecule = file_io.get_molecule_from_file_fchk(data_dir + '/file_io/test.fchk', read_multiple=False)
        self.assertEqual(molecule.electronic_structure.alpha_occupancy, alpha_occupancy)
        self.assertEqual(molecule.electronic_structure.beta_occupancy, beta_occupancy)
        self.assertEqual(molecule.electronic_structure.alpha_electrons, alpha_electrons)
        self.assertEqual(molecule.electronic_structure.beta_electrons, beta_electrons)

    def test_readinf_molden(self):
        alpha_occupancy = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        beta_occupancy = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        alpha_electrons = 13
        beta_electrons = 13
        molecule = file_io.get_molecule_from_file_molden(data_dir + '/file_io/test.molden', read_multiple=False)
        self.assertEqual(molecule.electronic_structure.alpha_occupancy, alpha_occupancy)
        self.assertEqual(molecule.electronic_structure.beta_occupancy, beta_occupancy)
        self.assertEqual(molecule.electronic_structure.alpha_electrons, alpha_electrons)
        self.assertEqual(molecule.electronic_structure.beta_electrons, beta_electrons)