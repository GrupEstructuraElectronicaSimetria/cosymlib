import numpy as np
from warnings import warn
from cosymlib.symmetry.pointgroup.operations import rotation_matrix


def print_symmetry_labels():
    operations = {'E': 'Identity Symmetry',
                  'Ci': 'Inversion Symmetry Group',
                  'Cs': 'Reflection Symmetry Group',
                  'Cn': 'Rotational Symmetry Group (n: rotation order)',
                  'Sn': 'Rotation-Reflection Symmetry Group (n: rotation-reflection order)'}

    print('Available symmetry groups\n')
    for k, label in operations.items():
        print('{:6} {}'.format(k, label))

# def get_group_num_from_label(label):
#
#     operations = {'Cn':  '',
#                   'CnH': [2, ],
#                   'CnV': [3, ],
#                   'Ci':  [0, 2],
#                   'Cs':  [0, 3],
#                   'Cinf':[9, 1],
#                   'Dn':  [4, ],
#                   'DnH': [5, ],
#                   'DnD': [6, ],
#                   'Dinf':[9, 2],
#                   'Sn':  [7, ],
#                   'T':   [8, 1],
#                   'Th':  [8, 2],
#                   'Td':  [8, 3],
#                   'O':   [8, 4],
#                   'Oh':  [8, 5],
#                   'I':   [8, 6],
#                   'Ih':  [8, 7],
#                   }


def orthogonal_c4(c3, c4):
    def calculate_angle(axis1, axis2):
        axis1 /= np.linalg.norm(axis1)
        axis2 /= np.linalg.norm(axis2)
        return np.arccos(np.dot(axis1, axis2) / (np.linalg.norm(axis1) * np.linalg.norm(axis2))) * 180 / np.pi

    perfect_c3_c4 = 125.26438968275465
    c3 /= np.linalg.norm(c3)
    c4 /= np.linalg.norm(c4)
    angle_c3_c4 = calculate_angle(c3, c4)
    diff_angle = perfect_c3_c4 - angle_c3_c4
    if abs(diff_angle) <= 1E-8:
        new_c4 = c4
    else:
        new_c4 = np.dot(rotation_matrix(np.cross(c3, c4), diff_angle), c4)
        warn('ChangedAxisWarning: Set C4 axis to: {}'.format(new_c4))

    rot_mat = rotation_matrix(c3, 120.0)
    c4_2 = np.dot(rot_mat, new_c4)
    return new_c4, c4_2/np.linalg.norm(c4_2)