from cosymlib.shape import tools, Shape
import numpy as np
from cosymlib.molecule.geometry import Geometry

def get_shape_map(shape_label1, shape_label2, num_points=20):

    path_structures = []
    if isinstance(shape_label1, Geometry):
        coordinates_a = shape_label1.get_positions()
    else:
        coordinates_a = tools.get_test_structure(shape_label1, central_atom=1).get_positions()
    coordinates_b = Shape(coordinates_a).structure(shape_label2, central_atom=len(coordinates_a))
    coordinates_b = np.concatenate((coordinates_b[1:], [coordinates_b[0]]))


    measure_label1 = []
    measure_label2 = []
    for i in range(0, num_points+1):
        structure_along_path = []
        for ida, atom in enumerate(coordinates_a):
            structure_along_path.append(atom + (coordinates_b[ida] - atom)*i/num_points)

        measure_label1.append(Shape(structure_along_path).measure(shape_label1, central_atom=len(coordinates_a)))
        measure_label2.append(Shape(structure_along_path).measure(shape_label2, central_atom=len(coordinates_a)))
        path_structures.append(structure_along_path)

    return measure_label1, measure_label2, path_structures
