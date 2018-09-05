import numpy as np
from symeess.molecule.geometry import Geometry
from symeess.shape import get_test_structure


def get_shape_map(shape_label1, shape_label2, central_atom=None, num_points=50):
    ideal_structure = get_test_structure(shape_label1, central_atom)
    ideal_label_structure = []
    for idx, atom in enumerate(ideal_structure):
        if idx == central_atom-1:
            atom = np.ndarray.tolist(atom)
            atom.insert(0, 'M')
        else:
            atom = np.ndarray.tolist(atom)
            atom.insert(0, 'L')
        ideal_label_structure.append(atom)

    geometry = Geometry(ideal_label_structure)
    S_label1 = [0]
    S_label2 = [geometry.get_shape_measure(shape_label2, central_atom=central_atom)]
    try:
        theta = minimum_distortion_angles[shape_label1][shape_label2]
    except KeyError:
        theta = minimum_distortion_angles[shape_label2][shape_label1]
    dtheta = np.linspace(0, theta, num_points)
    for angle in dtheta[1:]:
        a_label1 = angle
        a_label2 = theta - a_label1
        S_label1.append(100 * (np.sin(np.radians(a_label1)) ** 2))
        S_label2.append((100 * (np.sin(np.radians(a_label2)) ** 2)))
    return S_label1, S_label2


def get_path_deviation(Sx, Sy, shape_label1, shape_label2):
    path_deviation = []
    for idx, x in enumerate(Sx):
        new_theta = np.arcsin(np.sqrt(x)/10) + np.arcsin(np.sqrt(Sy[idx])/10)
        try:
            theta = minimum_distortion_angles[shape_label1][shape_label2]
        except KeyError:
            theta = minimum_distortion_angles[shape_label2][shape_label1]
        path_deviation.append(((new_theta/np.radians(theta))-1)*100)

    return path_deviation

minimum_distortion_angles = {'T-4': {'SS-4': 18.234, 'SP-4': 35.264}, 'SS-4': {'SP-4': 25.878},
                             'vOC-5': {'SPY-5': 7.582, 'TBPY-5': 15.722, 'PP-5': 34.588},
                             'SPY-5': {'TBPY-5': 13.417, 'PP-5': 35.243}, 'TBPY-5': {'PP-5': 37.506},
                             'OC-6': {'TPR-6': 24.149, 'PPY-6': 33.484, 'HP-6': 35.264},
                             'TPR-6': {'PPY-6': 24.362, 'HP-6': 35.472}, 'PPY-6': {'HP-6': 32.359},
                             'COC-7': {'CTPR-7': 7.099, 'JPBPY-7': 16.852, 'HPY-7': 24.393, 'HP-7': 37.924},
                             'CTPR-7': {'JPBPY-7': 14.934, 'HPY-7': 26.530, 'HP-7': 36.794},
                             'JPBPY-7': {'HPY-7': 31.105, 'HP-7': 36.399}, 'HPY-7': {'HP-7': 30.309},
                             'CU-8': {'TDD-8': 16.379, 'SAPR-8': 19.360, 'HBPY-8': 16.842, 'HPY-8': 33.592,
                                      'OP-8': 38.240}, 'TDD-8': {'SAPR-8': 9.716, 'HBPY-8': 23.326, 'HPY-8': 29.863,
                             'OP-8': 34.533}, 'SAPR-8': {'HBPY-8': 25.444, 'HPY-8': 29.691, 'OP-8': 30.736},
                             'HBPY-8': {'HPY-8': 29.109, 'OP-8': 34.708}, 'HPY-8': {'OP-8': 28.528}}
