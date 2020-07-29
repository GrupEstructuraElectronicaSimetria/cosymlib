import numpy as np
from cosymlib.file_io import read_generic_structure_file, get_file_xyz_txt
from cosymlib.molecule.geometry import Geometry


def read_normal_modes_gaussian(file_name):
    read = False
    frequencies = []
    normal_modes_matrices = []
    with open(file_name, mode='r') as lines:
        for line in lines:
            if read:
                if len(line.split()) == 3 or not line.strip():
                    read = False
                else:
                    for i in range(int(len(line.split()[2:])/3)):
                        atom_vibration = line.split()[2 + 3*i:3 + 2 + 3*i]
                        normal_modes_matrices[-3+i].append([float(x)for x in atom_vibration])

            if 'Atom  AN' in line:
                read = True
                for _ in range(line.split()[2:].count('X')):
                    normal_modes_matrices.append([])
            if 'Frequencies --' in line:
                frequencies.append(line.split()[2:])
    frequencies = [float(frequency) for sublist in frequencies for frequency in sublist]

    return frequencies, normal_modes_matrices


def get_nm_vibration_path(geometry, normal_mode, points=10, backward=False, k=1.5):

    nm_np_matrix = np.array(normal_mode)
    geom_np = np.array(geometry)
    if backward:
        x = np.linspace(-k, k, points)
    else:
        x = np.linspace(k/points, k, points)

    structures_path = []
    for dx in x:
        structures_path.append(geom_np + dx*nm_np_matrix)

    return structures_path


molecule = read_generic_structure_file('sf6.fchk')
freq, nm_martices = read_normal_modes_gaussian('sf6_freq.out')
k_points = 1.5
points = 10
n_freq = 1
path = get_nm_vibration_path(molecule.geometry.get_positions(), nm_martices[n_freq - 1], k=k_points)

x = np.linspace(k_points/points, k_points, points)
output = open('sf6_freq{}.xyz'.format(n_freq), 'w')
for ids, structure in enumerate(path):
    xi = '{:.2f}'.format(x[ids]).replace('.', '')
    txt = get_file_xyz_txt(Geometry(structure, symbols=molecule.geometry.get_symbols(), name='freq : ' + str(freq[0]) +
                                                                                             ' ' + str(xi)))
    output.write(txt)
