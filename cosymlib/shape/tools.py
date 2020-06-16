import os
import yaml
import numpy as np

ideal_structures = None


def get_test_structure(label, central_atom=0):
    global ideal_structures
    from cosymlib.molecule.geometry import Geometry

    if ideal_structures is None:
        file_path = os.path.dirname(os.path.abspath(__file__)) + '/ideal_structures_center.yaml'
        with open(file_path, 'r') as stream:
            ideal_structures = yaml.load(stream, Loader=yaml.FullLoader)

    if central_atom == 0:
        coordinates = ideal_structures[label][:-1]
        return Geometry(positions=coordinates,
                        name=label,
                        symbols=['L'] * len(coordinates),
                        connectivity=[])

    else:
        coordinates = ideal_structures[label]
        return Geometry(positions=coordinates,
                        name=label,
                        symbols=['L'] * (len(coordinates)-1) + ['M'],
                        connectivity=[])


def get_structure_references(vertices):

    # print('Available reference structure with {} Vertices'.format(vertices))
    references_list = []
    for ref in shape_structure_references['{} Vertices'.format(vertices)]:
        references_list.append(ref[0])
    return references_list


def order_coordinates(coordinates, indices):
    indices = [i - 1 for i in indices]
    coordinates = np.array(coordinates).copy()
    coordinates[indices, :] = coordinates[indices[::-1], :]
    return coordinates


def get_shape_label(code, vertices):

    for label in shape_structure_references[str(vertices)+' Vertices']:
        if label[1] == code:
            return label[0]


def get_shape_label_info(vertices, old=False, with_central_atom=False):
    if with_central_atom:
        vertices = vertices-1

    txt = 'Available reference structures with {} Vertices:\n\n'.format(vertices)
    if old:
        txt += '{:<4} '.format('')
    txt += '{:10}  {:8}  {}\n\n'.format('Label', 'Sym', 'Info')
    for idl, labels in enumerate(shape_structure_references['{} Vertices'.format(vertices)]):
        if old:
            txt += '{:<4} '.format(idl+1)
        txt += '{:10}  {:8}  {}\n'.format(labels[0], labels[2], labels[3])
    return txt


shape_structure_references = {'2 Vertices': [['L-2', 1, 'Dinfh', 'Linear'],
                                             ['vT-2', 2, 'C2v', 'Divacant tetrahedron'],
                                             ['vOC-2', 3, 'C2v', 'Tetravacant octahedron']],
                              '3 Vertices': [['TP-3', 1, 'D3h', 'Trigonal'],
                                             ['vT-3', 2, 'C3v', 'Vacant tetrahedron'],
                                             ['fvOC-3', 3, 'C3v', 'fac-Trivacant octahedron'],
                                             ['mvOC-3', 4, 'C2v', 'mer-Trivacant octahedron']],
                              '4 Vertices': [['SP-4', 1, 'D4h', 'Square'],
                                             ['T-4', 2, 'Td', 'Tetrahedron'],
                                             ['SS-4', 3, 'C2v', 'Seesaw'],
                                             ['vTBPY-4', 4, 'C3v', 'Axially vacant trigonal bipyramid']],
                              '5 Vertices': [['PP-5', 1, 'D5h', 'Pentagon'],
                                             ['vOC-5', 2, 'C4v', 'Vacant octahedron'],
                                             ['TBPY-5', 3, 'D3h', 'Trigonal bipyramid'],
                                             ['SPY-5', 4, 'C4v', 'Spherical square pyramid'],
                                             ['JTBPY-5', 5, 'D3h', 'Johnson trigonal bipyramid J12']],
                              '6 Vertices': [['HP-6', 1, 'D6h', 'Hexagon'],
                                             ['PPY-6', 2, 'C5v', 'Pentagonal pyramid'],
                                             ['OC-6', 3, 'Oh', 'Octahedron'],
                                             ['TPR-6', 4, 'D3h', 'Trigonal prism'],
                                             ['JPPY-6', 5, 'C5v', 'Johnson pentagonal pyramid J2']],
                              '7 Vertices': [['HP-7', 1, 'D7h', 'Heptagon'],
                                             ['HPY-7', 2, 'C6v', 'Hexagonal pyramid'],
                                             ['PBPY-7', 3, 'D5h', 'Pentagonal bipyramid'],
                                             ['COC-7', 4, 'C3v', 'Capped octahedron'],
                                             ['CTPR-7', 5, 'C2v', 'Capped trigonal prism'],
                                             ['JPBPY-7', 6, 'D5h', 'Johnson pentagonal bipyramid J13'],
                                             ['JETPY-7', 7, 'C3v', 'Johnson elongated triangular pyramid J7']],
                              '8 Vertices': [['OP-8', 1, 'D8h', 'Octagon'],
                                             ['HPY-8', 2, 'C7v', 'Heptagonal pyramid'],
                                             ['HBPY-8', 3, 'D6h', 'Hexagonal bipyramid'],
                                             ['CU-8', 4, 'Oh', 'Cube'],
                                             ['SAPR-8', 5, 'D4d', 'Square antiprism'],
                                             ['TDD-8', 6, 'D2d', 'Triangular dodecahedron'],
                                             ['JGBF-8', 7, 'D2d', 'Johnson gyrobifastigium J26'],
                                             ['JETBPY-8', 8, 'D3h', 'Johnson elongated triangular bipyramid J14'],
                                             ['JBTPR-8', 9, 'C2v', 'Biaugmented trigonal prism J50'],
                                             ['BTPR-8', 10, 'C2v', 'Biaugmented trigonal prism'],
                                             ['JSD-8', 11, 'D2d', 'Snub diphenoid J84'],
                                             ['TT-8', 12, 'Td', 'Triakis tetrahedron'],
                                             ['ETBPY-8', 13, 'D3h', 'Elongated trigonal bipyramid']],
                              '9 Vertices': [['EP-9', 1, 'D9h', 'Enneagon'],
                                             ['OPY-9', 2, 'C8v', 'Octagonal pyramid'],
                                             ['HBPY-9', 3, 'D7h', 'Heptagonal bipyramid'],
                                             ['JTC-9', 4, 'C3v', 'Johnson triangular cupola J3'],
                                             ['JCCU-9', 5, 'C4v', 'Capped cube J8'],
                                             ['CCU-9', 6, 'C4v', 'Spherical-relaxed capped cube'],
                                             ['JCSAPR-9', 7, 'C4v', 'Capped square antiprism J10'],
                                             ['CSAPR-9', 8, 'C4v', 'Spherical capped square antiprism'],
                                             ['JTCTPR-9', 9, 'D3h', 'Tricapped trigonal prism J51'],
                                             ['TCTPR-9', 10, 'D3h', 'Spherical tricapped trigonal prism'],
                                             ['JTDIC-9', 11, 'C3v', 'Tridiminished icosahedron J63'],
                                             ['HH-9', 12, 'C2v', 'Hula-hoop'],
                                             ['MFF-9', 13, 'Cs', 'Muffin']],
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
                              '48 Vertices': [['TCOC-48', 1, 'Oh', 'Truncated cuboctahedron']],
                              '60 Vertices': [['TRIC-60', 1, 'Ih', 'Truncated icosahedron (fullerene)']]}