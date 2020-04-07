def print_symmetry_labels():
    operations = {'E': 'Identity Symmetry',
                  'Ci': 'Inversion Symmetry Group',
                  'Cs': 'Reflection Symmetry Group',
                  'Cn': 'Rotational Symmetry Group (n: rotation order)',
                  'Sn': 'Rotation-Reflection Symmetry Group (n: rotation-reflection order)'}

    print('Available symmetry groups\n')
    for k, label in operations.items():
        print('{:6} {}'.format(k, label))

def get_group_num_from_label(label):

    operations = {'Cn':  '',
                  'CnH': [2, ],
                  'CnV': [3, ],
                  'Ci':  [0, 2],
                  'Cs':  [0, 3],
                  'Cinf':[9, 1],
                  'Dn':  [4, ],
                  'DnH': [5, ],
                  'DnD': [6, ],
                  'Dinf':[9, 2],
                  'Sn':  [7, ],
                  'T':   [8, 1],
                  'Th':  [8, 2],
                  'Td':  [8, 3],
                  'O':   [8, 4],
                  'Oh':  [8, 5],
                  'I':   [8, 6],
                  'Ih':  [8, 7],
                  }
