
def write_shape_measure_data(measures, molecules_name, shape_label):

    txt_shape = '{}'.format('Structure')
    max_name = len(max(molecules_name, key=len))
    if max_name < 9:
        n = 5
    else:
        n = max_name - 4
    for label in shape_label:
        n += len(label)
        txt_shape += '{}'.format(label.rjust(n))
        n = 12 - len(label)
    txt_shape += '\n\n'

    for idx, name in enumerate(molecules_name):
        max_name = len(max(molecules_name, key=len))
        txt_shape += '{}'.format(name)
        if max_name < 9:
            n = 18 - len(name)
        else:
            n = 9 + max_name - len(name)
        for idn, label in enumerate(shape_label):
            txt_shape += ' {:{width}.{prec}f} ,'.format(measures[idn][idx], width=n, prec=3)
            n = 11
        txt_shape += '\n'
    txt_shape += '\n'

    return txt_shape


# def write_shape_structure_data(geometries, structures, measures, shape_label, output_name=None):
#
#     if output_name is not None:
#         output = open(output_name + '.out', 'w')
#     else:
#         output = sys.stdout
#     shape_header(output)
#
#     for idx, geometry in enumerate(geometries):
#         output.write('\n')
#         output.write('Structure {} : {}\n'.format(idx+1, geometry.get_name()))
#
#         for idn, array in enumerate(geometry.get_positions()):
#             output.write('{:2s}'.format(geometry.get_symbols()[idn]))
#             output.write(' {:11.7f} {:11.7f} {:11.7f}\n'.format(array[0], array[1], array[2]))
#         output.write('\n')
#
#         for idn, label in enumerate(shape_label):
#             output.write('{} Ideal Structure CShM = {:.3f}\n'
#                          .format(label, measures[idn][idx]))
#             for jd, array in enumerate(structures[idn][idx]):
#                 output.write('{:2s}'.format(geometry.get_symbols()[jd]))
#                 output.write(' {:11.7f} {:11.7f} {:11.7f}\n'.format(array[0], array[1], array[2]))
#             output.write('\n')
#
#         output.write('-' * 70 + '\n')
    # output.close()


def write_shape_map(shape_label1, shape_label2, path):

    txt_shape = " {:6} {:6}\n".format(shape_label1, shape_label2)
    for idx, value in enumerate(path[0]):
        txt_shape += '{:6.3f}, {:6.3f}'.format(path[0][idx], path[1][idx])
        txt_shape += '\n'
    return txt_shape


def write_minimal_distortion_path_analysis(shape_label1, shape_label2, measures, pathdev, GenCoord,
                                           maxdev, mindev, mingco, maxgco, molecule_names):

    txt_shape = "Deviation threshold to calculate Path deviation function: {:2.1f}% - {:2.1f}%\n".format(mindev, maxdev)
    txt_shape += "Deviation threshold to calculate Generalized Coordinate: {:2.1f}% - {:2.1f}%\n".format(mingco, maxgco)
    txt_shape += "\n"
    txt_shape += '{:}'.format('structure'.upper())
    txt_shape += " {:>7} {:>9}".format(shape_label1, shape_label2)
    txt_shape += "{:>12} {:>9}".format('DevPath', 'GenCoord')
    txt_shape += "\n"

    for idx, molecule_name in enumerate(molecule_names):
        txt_shape += '{}  ,'.format(molecule_name)
        if molecule_names[idx].strip() == '':
            width = 6 + len(molecule_names[idx])
        else:
            width = 14 - len(molecule_names[idx])
        for label in [shape_label1, shape_label2]:
            txt_shape += ' {:{width}.{prec}f},'.format(measures[label][idx], width=width, prec=3)
            width = 7
        txt_shape += '{:8.1f}, {:8}'.format(pathdev[idx], GenCoord[idx])
        txt_shape += '\n'
    return txt_shape

