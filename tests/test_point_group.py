import unittest
from cosymlib import file_io, Cosymlib
import os


data_dir = os.path.join(os.path.dirname(__file__), 'data')


class TestPointGroup(unittest.TestCase):

    def test(self):
        point_groups = ['C1', 'Cs', 'Ci', 'Cinfh', 'Dinfh', 'C2', 'C3', 'C2h', 'C3h', 'C2v', 'C3v', 'C4v', 'C5v', 'D2',
                        'D3', 'D2h', 'D3h', 'D4h', 'D5h', 'D6h', 'D7h', 'D8h', 'D2d', 'D3d', 'D4d', 'D5d', 'D6d', 'D8d',
                        'S4', 'T', 'Th', 'Td', 'O', 'Oh']
        molecules = file_io.read_generic_structure_file(data_dir + '/point_group/sym_molecules.xyz', read_multiple=True)
        symobj = Cosymlib(molecules)
        calculated_point_groups = [pg for pg in (symobj.get_point_group())]
        for pg1, pg2 in zip(point_groups, calculated_point_groups):
            self.assertEqual(pg1, pg2)
