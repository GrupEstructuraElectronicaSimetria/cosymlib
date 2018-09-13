import os
import sys


def write_shape_data(data, shape_label, molecule_names, option, output_name=sys.stdout):

    extensions = {'measure': '.tab', 'structure': '.out', 'test': '.tst'}
    if not os.path.exists('./results'):
        os.makedirs('./results')
    output = open('results/'+output_name + extensions[option], 'w')

    output.write('-'*40 + '\n')
    output.write('Shape measure/s \n')
    output.write('-'*40 + '\n')

    if 'measure' in option:
        output.write('{}'.format('structure'.upper()))
        for label in shape_label:
            n = len(label) + 3
            output.write('{}'.format(label.rjust(n)))
        output.write('\n')
        for idx, molecule_name in enumerate(molecule_names):
            output.write('{}'.format(molecule_name))
            if molecule_names[idx].strip() == '':
                n = 4 + len(molecule_names[idx])
            else:
                n = 14 - len(molecule_names[idx])
            for label in shape_label:
                output.write(' {:{width}.{prec}f}'.format(data[label][idx], width=n, prec=3))
                n = 7
            output.write('\n')
        output.write('\n')

    if 'structure' in option:
        output.write("{}\n".format('ideal_structure'.upper()))
        for idx, molecule_name in enumerate(molecule_names):
            output.write('\n')
            output.write('{}'.format(molecule_names[idx]))

            n = int(23 - len(molecule_names[idx]))
            for label in shape_label:
                output.write('{}'.format(label.rjust(n)))
                n = 37
            output.write('\n')

            for idn, symbol in enumerate(data['symbols'][idx]):
                output.write('{:2s}'.format(symbol))
                for label in shape_label:
                    array = data[label][idx][idn]
                    output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
                output.write('\n')

        output.write('\n')

    if 'test_structure' in option:
        output.write("{}\n".format('test_structure'.upper()))
        n = 20
        for label in shape_label:
            output.write('{}'.format(label.rjust(n)))
            n = 36 + len(label)
        output.write('\n')

        for idx in list(range(len(data[0]['symbols']))):
            for label in shape_label:
                array = data[0][label]['test_structure'][idx]
                output.write(' {:11.7f} {:11.7f} {:11.7f} |'.format(array[0], array[1], array[2]))
            output.write('\n')

    output.close()


def write_shape_map(shape_label1, shape_label2, path, output_name='symeess_shape_map'):
    output = open('results/' + output_name + '.pth', 'w')
    output.write(" {:6} {:6}\n".format(shape_label1, shape_label2))
    for idx, value in enumerate(path[0]):
        output.write('{:6.3f} {:6.3f}'.format(path[0][idx], path[1][idx]))
        output.write('\n')


def write_minimal_distortion_path_analysis(shape_label1, shape_label2, measures, pathdev, GenCoord,
                                           maxdev, mindev, mingco, maxgco, molecule_names,
                                           output_name='symeess_shape'):
    output = open('results/' + output_name + '.flt', 'w')
    output.write("Deviation threshold to calculate Path deviation function: "
                 "{:2.1f}% - {:2.1f}%\n".format(mindev, maxdev))
    output.write("Deviation threshold to calculate Generalized Coordinate: "
                 "{:2.1f}% - {:2.1f}%\n".format(mingco, maxgco))
    output.write("\n")
    output.write('{:11}'.format('structure'.upper()))
    output.write(" {:7} {:7}".format(shape_label1, shape_label2))
    output.write("{:6} {:6}".format('DevPath', 'GenCoord'))
    output.write("\n")

    for idx, molecule_name in enumerate(molecule_names):
        output.write('{}'.format(molecule_name))
        if molecule_names[idx].strip() == '':
            n = 4 + len(molecule_names[idx])
        else:
            n = 15 - len(molecule_names[idx])
        for label in [shape_label1, shape_label2]:
            output.write(' {:{width}.{prec}f}'.format(measures[label][idx], width=n, prec=3))
            n = 7
        output.write('{:8.1f} {:8}'.format(pathdev[idx], GenCoord[idx]))
        output.write('\n')