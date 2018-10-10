#!/usr/bin/env python
import argparse
import symeess
import sys


def main(argv):

    parser = argparse.ArgumentParser(description='Symeess')
    parser.add_argument('-input_file', type=str, default=None, help='input file name(+extension)')
    parser.add_argument('-o', '-output', dest='output_name', default=None, help='customize output file name')

    # Shape input flags
    group_shape = parser.add_argument_group('Shape', 'Shape_options')
    group_shape.add_argument('-m', '--measure',
                             action='store_true',
                             default=False,
                             help='return shape measure of input structure with reference polyhedra')
    group_shape.add_argument('-s', '--structure',
                             action='store_true',
                             default=False,
                             help='calculate the ideal structure for the input structure')
    group_shape.add_argument('-t', '--test',
                             action='store_true',
                             default=False,
                             help='print the reference structure of the given label')
    group_shape.add_argument('-c', action='store',
                             type=int,
                             default=None,
                             help='position of the central atom if exist')
    group_shape.add_argument('-label',
                             dest='reference_polyhedra',
                             action='store',
                             nargs='*',
                             default=None,
                             help='use labels from Shape manual for desire reference polyhedra')
    group_shape.add_argument('-n',
                             action='store',
                             type=str,
                             default=None,
                             help='Print all the possible reference structures of n vertices')
    group_shape.add_argument('-map',
                             action='store_true',
                             default=False)
    group_shape.add_argument('-path',
                             action='store_true',
                             default=False)

    args = parser.parse_args(argv)

    # Reading and initializing
    if args.reference_polyhedra is not None:
        reference_polyhedron = args.reference_polyhedra # .split()
    if args.input_file is not None:
        geometries = symeess.file_io.read_input_file(args.input_file)
        example = symeess.Symeess()
        example.set_molecules(geometries)
    if args.output_name is not None:
        output_name = args.output_name
    else:
        output_name = 'symeess'
    central_atom = args.c

    # Shape's commands
    if args.structure:
        example.write_shape_structure_2file(reference_polyhedron, central_atom=central_atom, output_name=output_name)
    if args.measure:
        example.write_shape_measure_2file(reference_polyhedron, central_atom=central_atom, output_name=output_name)
    if args.test:
        for reference in reference_polyhedron:
            print(reference)
            for array in symeess.shape.shape_tools.get_test_structure(reference, central_atom=central_atom):
                print('{:11.8f} {:11.8f} {:11.8f}'.format(array[0], array[1], array[2]))
    if args.n:
        for label in symeess.shape.shape_tools.get_structure_references(args.n):
            print(label)
    if args.map:
        symeess.write_minimum_distortion_path_shape_2file(reference_polyhedron[0],
                                                          reference_polyhedron[1],
                                                          num_points=20)
    if args.path:
        example.write_path_parameters_2file('SP-4', 'T-4', central_atom=central_atom, output_name=output_name)

    # input("Press enter to exit")


if __name__ == '__main__':

    if sys.argv[1:]:
        argv = sys.argv[1:]
    else:
        argv = (['-m', '-label', 'SP-4', 'T-4',
                 '-s',
                 '-c', '1',
                 '-o', 'coord',
                 '-input_file', '../examples/coord.xyz'])
    main(argv)
