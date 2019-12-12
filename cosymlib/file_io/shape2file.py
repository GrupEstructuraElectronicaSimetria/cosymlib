
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
            txt_shape += ', {:{width}.{prec}f}'.format(measures[idn][idx], width=n, prec=3)
            n = 11
        txt_shape += '\n'
    txt_shape += '\n'

    return txt_shape


def write_shape_map(shape_label1, shape_label2, path):

    txt_shape = " {:6} {:6}\n".format(shape_label1, shape_label2)
    for idx, value in enumerate(path[0]):
        txt_shape += '{:6.3f}, {:6.3f}'.format(path[0][idx], path[1][idx])
        txt_shape += '\n'
    return txt_shape


def write_minimal_distortion_path_analysis(measures, pathdev, GenCoord, maxdev, mindev, mingco, maxgco, molecules_name):

    txt_shape = "Deviation threshold to calculate Path deviation function: {:2.1f}% - {:2.1f}%\n".format(mindev, maxdev)
    txt_shape += "Deviation threshold to calculate Generalized Coordinate: {:2.1f}% - {:2.1f}%\n".format(mingco, maxgco)
    txt_shape += "\n"
    txt_shape += '{:}'.format('structure'.upper())
    txt_shape += " {:>7} {:>9}".format(list(measures.keys())[0], list(measures.keys())[1])
    txt_shape += "{:>12} {:>9}".format('DevPath', 'GenCoord')
    txt_shape += "\n"

    for idx, molecule_name in enumerate(molecules_name):
        txt_shape += '{}  ,'.format(molecule_name)
        if molecule_name.strip() == '':
            width = 6 + len(molecule_name)
        else:
            width = 14 - len(molecule_name)
        for label in list(measures.keys()):
            txt_shape += ' {:{width}.{prec}f},'.format(measures[label][idx], width=width, prec=3)
            width = 7
        txt_shape += '{:8.1f}, {:8}'.format(pathdev[idx], GenCoord[idx])
        txt_shape += '\n'
    return txt_shape

