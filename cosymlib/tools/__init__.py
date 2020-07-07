import os
import yaml
import numpy as np
from warnings import warn


file_path = os.path.dirname(os.path.abspath(__file__)) + '/periodic_table.yaml'

with open(file_path, 'r') as stream:
    periodic_table_info = yaml.safe_load(stream)


def atomic_number_to_element(z):
    for element, data in periodic_table_info.items():
        if int(z) == data['atomic_number']:
            return element


def element_to_atomic_number(symbol):
    return periodic_table_info[symbol.capitalize()]['atomic_number']


def element_mass(symbol):
    return periodic_table_info[symbol.capitalize()]['atomic_mass']


def element_valence_electron(symbol):
    return periodic_table_info[symbol.capitalize()]['valence_electrons']


def center_mass(elements, coordinates):
    cm = [0., 0., 0.]
    m = 0
    for ide, element in enumerate(elements):
        cm += coordinates[ide] * element_mass(element)
        m += element_mass(element)
    cm = cm/m
    return cm


def generate_connectivity_from_geometry_slow(geometry, thresh=1.2):
    coor = geometry.get_positions()
    sym = geometry.get_symbols()

    connectivity = []
    for i, (atom1, sym1) in enumerate(zip(np.array(coor), sym)):
        for j, (atom2, sym2) in enumerate(zip(np.array(coor), sym)):
            dist = np.linalg.norm(atom1 - atom2)
            rad = periodic_table_info[sym1]['covalent_radius'] + periodic_table_info[sym2]['covalent_radius']

            if np.abs(rad - dist)/rad < thresh - 1:
                connectivity.append((i+1, j+1))

    return connectivity


def generate_connectivity_from_geometry(geometry, thresh=1.2):
    from scipy.spatial import distance_matrix
    coordinates = geometry.get_positions()
    try:
        radius = [periodic_table_info[sym]['covalent_radius'] for sym in geometry.get_symbols()]
    except KeyError:
        warn('failed to generate connectivity, no connectivity will be used')
        return None

    distances_matrix = distance_matrix(coordinates, coordinates)

    radii_matrix = np.array([radius]*len(radius))
    radii_matrix = radii_matrix + radii_matrix.T

    try:
        relative_differences = np.abs(radii_matrix - distances_matrix)/radii_matrix
    except ValueError:
        warn('failed to generate connectivity')
        return None

    return (np.array(np.where(relative_differences < thresh - 1)).T + 1).tolist()


def get_connectivity_matrix(connectivity, ndim):
    cm = np.zeros((ndim, ndim), dtype=int)
    for pair in connectivity:
        cm[pair[0]-1, pair[1]-1] = 1
    return cm
