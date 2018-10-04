# import numpy as np
# from symeess.molecule.geometry import Geometry
# from symeess.shape.shape_tools import get_test_structure, _get_symmetry_angle
from symeess.shape import shape_tools, Shape


# def get_shape_map_old(shape_label1, shape_label2, central_atom=None, num_points=50):
#     ideal_structure = get_test_structure(shape_label1)
#     ideal_label_structure = []
#     symbols = []
#     for idx, atom in enumerate(ideal_structure):
#         if idx == 0 and central_atom is not None:
#             atom = np.ndarray.tolist(atom)
#             symbols.append('M')
#         else:
#             atom = np.ndarray.tolist(atom)
#             symbols.append('L')
#         ideal_label_structure.append(atom)
#
#     geometry = Geometry(symbols=symbols,
#                         positions=ideal_label_structure)
#     S_label1 = [0]
#     S_label2 = [geometry.get_shape_measure(shape_label2)]
#     theta = _get_symmetry_angle(shape_label1, shape_label2)
#     dtheta = np.linspace(0, theta, num_points)
#     for angle in dtheta[1:]:
#         a_label1 = angle
#         a_label2 = theta - a_label1
#         S_label1.append(100 * (np.sin(np.radians(a_label1)) ** 2))
#         S_label2.append((100 * (np.sin(np.radians(a_label2)) ** 2)))
#     return S_label1, S_label2


def get_shape_map(shape_label1, shape_label2, central_atom=None, num_points=20):

    structure_a = shape_tools.get_test_structure(shape_label1, central_atom=central_atom)
    structure_b = Shape(structure_a).structure(shape_label2, central_atom=5)
    structure_b = shape_tools.order_coordinates(structure_b, [1, len(structure_b)])

    measure_label1 = []
    measure_label2 = []
    for i in range(0, num_points+1):
        structure_along_path = []
        for ida, atom in enumerate(structure_a):
            structure_along_path.append(atom + (structure_b[ida] - atom)*i/num_points)

        measure_label1.append(Shape(structure_along_path).measure(shape_label1, central_atom=5))
        measure_label2.append(Shape(structure_along_path).measure(shape_label2, central_atom=5))

    # print(shape_label1, shape_label2)
    # for idm, measure in enumerate(measure_label1):
    #     print('{:.3f}  {:.3f}'.format(measure, measure_label2[idm]))
    return measure_label1, measure_label2
