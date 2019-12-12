from cosymlib.shape import shape_tools, Shape
import numpy as np


def get_shape_map(shape_label1, shape_label2, num_points=20):

    path_structures = []
    if isinstance(shape_label1, np.ndarray):
        structure_a = shape_label1
    else:
        structure_a = shape_tools.get_test_structure(shape_label1, central_atom=1)
    structure_b = Shape(structure_a).structure(shape_label2, central_atom=len(structure_a))
    structure_b = np.concatenate((structure_b[1:], [structure_b[0]]))

    measure_label1 = []
    measure_label2 = []
    for i in range(0, num_points+1):
        structure_along_path = []
        for ida, atom in enumerate(structure_a):
            structure_along_path.append(atom + (structure_b[ida] - atom)*i/num_points)

        measure_label1.append(Shape(structure_along_path).measure(shape_label1, central_atom=len(structure_a)))
        measure_label2.append(Shape(structure_along_path).measure(shape_label2, central_atom=len(structure_a)))
        path_structures.append(structure_along_path)

    return measure_label1, measure_label2, path_structures
