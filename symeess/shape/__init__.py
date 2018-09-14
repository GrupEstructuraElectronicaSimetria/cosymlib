from symeess.shape import shp
import yaml
import numpy as np
import os


def get_measure(geometry, shape_label, central_atom=None):
    """
    Compute the shape measure of the given geometry

    :param geometry: this object contains information about the positions of the molecule as well as the information
                     needed to carry out the shape's module
    :return: difference between user's structure and the one's that is compared. While 0 is no difference and 100
             is completly different
    """
    if central_atom is not None:
        coordinates = _order_coordinates(geometry.get_positions(), central_atom)
        code = get_ideal_structure(shape_label, geometry.get_n_atoms() - 1)
        c_atom = True
    else:
        coordinates = geometry.get_positions()
        code = get_ideal_structure(shape_label, geometry.get_n_atoms())
        c_atom = False
    measure_number = shp.cshm(coordinates, code, c_atom)
    return measure_number


def get_structure(geometry, shape_label, central_atom=None):
    """
    Calculate the ideal structure of the given geometry from a reference structure

    :param geometry: same as before
    :return: ideal structure if user's structure had the compared structure's shape
    """
    if central_atom is not None:
        coordinates = _order_coordinates(geometry.get_positions(), central_atom)
        code = get_ideal_structure(shape_label, geometry.get_n_atoms() - 1)
        c_atom = True
    else:
        coordinates = geometry.get_positions()
        code = get_ideal_structure(shape_label, geometry.get_n_atoms())
        c_atom = False
    measure_structure = shp.poly(coordinates, code, c_atom)
    return measure_structure[1], measure_structure[0]


def get_test_structure(shape_label, central_atom):
    file_path = os.path.dirname(os.path.abspath(__file__))+'/ideal_structures.yaml'
    with open(file_path, 'r') as stream:
        ideal_structures = yaml.load(stream)
    if central_atom is None:
        ideal_structures[shape_label].pop(0)
    measure_structure = np.array(ideal_structures[shape_label])
    return measure_structure


def get_ideal_structure(symbol, n_atoms):
    n_vertices = str(n_atoms)+' Vertices'
    for structure in shape_structure_references[n_vertices]:
        if structure[0] == symbol:
            return structure[1]
    raise NameError('Wrong ideal structure. N vertices != N atoms')


def get_path_deviation(Sx, Sy, shape_label1, shape_label2):
    new_theta = np.arcsin(np.sqrt(Sx)/10) + np.arcsin(np.sqrt(Sy)/10)
    theta = _get_symmetry_angle(shape_label1, shape_label2)
    path_deviation = ((new_theta/np.radians(theta))-1)*100
    return path_deviation


def get_generalized_coordinate(Sq, shape_label1, shape_label2):
    theta = _get_symmetry_angle(shape_label1, shape_label2)
    GenCoord = round(100*np.arcsin(np.sqrt(Sq)/10)/np.radians(theta), 1)
    return GenCoord


def _order_coordinates(coordinates, c_atom):
        c_atom = c_atom - 1
        new_coordinates = []
        for idx, array in enumerate(coordinates):
            if idx == c_atom:
                new_coordinates.insert(0, array)
            else:
                new_coordinates.append(array)
        return new_coordinates


def _get_shape_references(number_vertices):
    return shape_structure_references[number_vertices+' Vertices']


def _get_symmetry_angle(shape_label1, shape_label2):
    try:
        theta = minimum_distortion_angles[shape_label1][shape_label2]
    except KeyError:
        theta = minimum_distortion_angles[shape_label2][shape_label1]
    return theta


def get_shape_label(code, n_atoms):
    vertices = str(n_atoms)+' Vertices'
    label = shape_structure_references[vertices][code-1]
    return label[0]


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

shape_structure_references = {'2 Vertices': [['L-2', 1, 'Dinfh', 'Linear'],
                                             ['vT-2', 2, 'C2v', 'Divacant tetrahedron'],
                                             ['vOC-2', 3, 'C2v', 'Tetravacant octahedron']],
                              '3 Vertices': [['TP-3', 1, 'D3h', 'Trigonal'], ['vT-3', 2, 'C3v', 'Vacant tetrahedron'],
                                             ['fvOC-3', 3, 'C3v', 'fac-Trivacant octahedron'],
                                             ['mvOC-3', 4, 'C2v', 'mer-Trivacant octahedron']],
                              '4 Vertices': [['SP-4', 1, 'D4h', 'Square'], ['T-4', 2, 'Td', 'Tetrahedron'],
                                             ['SS-4', 3, 'C2v', 'Seesaw']],
                              '5 Vertices': [['PP-5', 1, 'D5h', 'Pentagon'], ['vOC-5', 2, 'C4v', 'Vacant octahedron'],
                                             ['TBPY-5', 3, 'D3h', 'Trigonal bipyramid'],
                                             ['SPY-5', 4, 'C4v', 'Spherical square pyramid'],
                                             ['JTBPY-5', 5, 'D3h', 'Johnson trigonal bipyramid J12']],
                              '6 Vertices': [['HP-6', 1, 'D6h', 'Hexagon'], ['PPY-6', 2, 'C5v', 'Pentagonal pyramid'],
                                             ['OC-6', 3, 'Oh', 'Octahedron'], ['TPR-6', 4, 'D3h', 'Trigonal prism'],
                                             ['JPPY-6', 5, 'C5v', 'Johnson pentagonal pyramid J2']],
                              '7 Vertices': [['HP-7', 1, 'D7h', 'Heptagon'], ['HPY-7', 2, 'C6v', 'Hexagonal pyramid'],
                                             ['PBPY-7', 3, 'D5h', 'Pentagonal bipyramid'],
                                             ['COC-7', 4, 'C3v', 'Capped octahedron'],
                                             ['CTPR-7', 5, 'C2v', 'Capped trigonal prism'],
                                             ['JPBPY-7', 6, 'D5h', 'Johnson pentagonal bipyramid J13'],
                                             ['JETPY-7', 7, 'C3v', 'Johnson elongated triangular pyramid J7']],
                              '8 Vertices': [['OP-8', 1, 'D8h', 'Octagon'], ['HPY-8', 2, 'C7v', 'Heptagonal pyramid'],
                                             ['HBPY-8', 3, 'D6h', 'Hexagonal bipyramid'], ['CU-8', 4, 'Oh', 'Cube'],
                                             ['SAPR-8', 5, 'D4d', 'Square antiprism'],
                                             ['TDD-8', 6, 'D2d', 'Triangular dodecahedron'],
                                             ['JGBF-8', 7, 'D2d', 'Johnson gyrobifastigium J26'],
                                             ['JETBPY-8', 8, 'D3h', 'Johnson elongated triangular bipyramid J14'],
                                             ['JBTPR-8', 9, 'C2v', 'Biaugmented trigonal prism J50'],
                                             ['BTPR-8', 10, 'C2v', 'Biaugmented trigonal prism'],
                                             ['JSD-8', 11, 'D2d', 'Snub diphenoid J84'],
                                             ['TT-8', 12, 'Td', 'Triakis tetrahedron'],
                                             ['ETBPY-8', 13, 'D3h', 'Elongated trigonal bipyramid']],
                              '9 Vertices': [['EP-9', 1, 'D9h', 'Enneagon'], ['OPY-9', 2, 'C8v', 'Octagonal pyramid'],
                                             ['HBPY-9', 3, 'D7h', 'Heptagonal bipyramid'],
                                             ['JTC-9', 4, 'C3v', 'Johnson triangular cupola J3'],
                                             ['JCCU-9', 5, 'C4v', 'Capped cube J8'],
                                             ['CCU-9', 6, 'C4v', 'Spherical-relaxed capped cube'],
                                             ['JCSAPR-9', 7, 'C4v', 'Capped square antiprism J10'],
                                             ['CSAPR-9', 8, 'C4v', 'Spherical capped square antiprism'],
                                             ['JTCTPR-9', 9, 'D3h', 'Tricapped trigonal prism J51'],
                                             ['TCTPR-9', 10, 'D3h', 'Spherical tricapped trigonal prism'],
                                             ['JTDIC-9', 11, 'C3v', 'Tridiminished icosahedron J63'],
                                             ['HH-9', 12, 'C2v', 'Hula-hoop'], ['MFF-9', 13, 'Cs', 'Muffin']],
                              '10 Vertices': [['DP-10', 1, 'D10h', 'Decagon'],
                                              ['EPY-10', 2, 'C9v', 'Enneagonal pyramid'],
                                              ['OBPY-10', 3, 'D8h', 'Octagonal bipyramid'],
                                              ['PPR-10', 4, 'D5h', 'Pentagonal prism'],
                                              ['PAPR-10', 5, 'D5d', 'Pentagonal antiprism'],
                                              ['JBCCU-10', 6, 'D4h', 'Bicapped cube J15'],
                                              ['JBCSAPR-10', 7, 'D4d', 'Bicapped square antiprism J17'],
                                              ['JMBIC-10', 8, 'C2v', 'Metabidiminished icosahedron J62'],
                                              ['JATDI-10', 9, 'C3v', 'Augmented tridiminished icosahedron J64'],
                                              ['JSPC-10', 10, 'C2v', 'Sphenocorona J87'],
                                              ['SDD-10', 11, 'D2', 'Staggered Dodecahedron (2:6:2)'],
                                              ['TD-10', 12, 'C2v', 'Tetradecahedron (2:6:2)'],
                                              ['HD-10', 13, 'D4h', 'Hexadecahedron (2:6:2) or (1:4:4:1)']],
                              '11 Vertices': [['HP-11', 1, 'D11h', 'Hendecagon'],
                                              ['DPY-11', 2, 'C10v', 'Decagonal pyramid'],
                                              ['EBPY-11', 3, 'D9h', 'Enneagonal bipyramid'],
                                              ['JCPPR-11', 4, 'C5v', 'Capped pentagonal prism J9'],
                                              ['JCPAPR-11', 5, 'C5v', 'Capped pentagonal antiprism J11'],
                                              ['JAPPR-11', 6, 'C2v', 'Augmented pentagonal prism J52'],
                                              ['JASPC-11', 7, 'Cs', 'Augmented sphenocorona J87']],
                              '12 Vertices': [['DP-12', 1, 'D12h', 'Dodecagon'],
                                              ['HPY-12', 2, 'C11v', 'Hendecagonal pyramid'],
                                              ['DBPY-12', 3, 'D10h', 'Decagonal bipyramid'],
                                              ['HPR-12', 4, 'D6h', 'Hexagonal prism'],
                                              ['HAPR-12', 5, 'D6d', 'Hexagonal antiprism'],
                                              ['TT-12', 6, 'Td', 'Truncated tetrahedron'],
                                              ['COC-12', 7, 'Oh', 'Cuboctahedron'],
                                              ['ACOC-12', 8, 'D3h', 'Anticuboctahedron J27'],
                                              ['IC-12', 9, 'Ih', 'Icosahedron'],
                                              ['JSC-12', 10, 'C4v', 'Johnson square cupola J4'],
                                              ['JEPBPY-12', 11, 'D6h', 'Johnson elongated pentagonal bipyramid J16'],
                                              ['JBAPPR-12', 12, 'C2v', 'Biaugmented pentagonal prism J53'],
                                              ['JSPMC-12', 13, 'Cs', 'Sphenomegacorona J88']],
                              '20 Vertices': [['DD-20', 1, 'Ih', 'Dodecahedron']],
                              '24 Vertices': [['TCU-24', 1, 'Oh', 'Truncated cube'],
                                              ['TOC-24', 2, 'Oh', 'Truncated octahedron']],
                              '48 Vertices' : ['TCOC-48', 1, 'Oh', 'Truncated cuboctahedron'],
                              '60 Vertices': ['TRIC-60', 1, 'Ih', 'Truncated icosahedron (fullerene)']}
