from cosymlib.shape import maps
import numpy as np
import sys


def plot_minimum_distortion_path_shape(shape_label1, shape_label2, num_points=20, output=sys.stdout, show_plot=True):
    import matplotlib.pyplot as plt

    path = get_shape_path(shape_label1, shape_label2, num_points)

    shape_map_txt = " {:6} {:6}\n".format(shape_label1, shape_label2)
    for idx, value in enumerate(path[0]):
        shape_map_txt += '{:6.3f}, {:6.3f}'.format(path[0][idx], path[1][idx])
        shape_map_txt += '\n'

    print(shape_map_txt)
    if show_plot:
        plt.plot(path[0], path[1], 'k', linewidth=2.0)
        plt.xlabel(shape_label1)
        plt.ylabel(shape_label2)
        plt.show()


def get_shape_path(shape_label1, shape_label2, num_points):
    return maps.get_shape_map(shape_label1, shape_label2, num_points)


def plot_molecular_orbital_diagram(molecule, wfnsym, mo_range=None):
    import matplotlib.pyplot as plt

    labels = wfnsym.IRLab
    if mo_range is not None:
        ird_a_max = [np.argmax(ird_a_orb) for ird_a_orb in wfnsym.mo_IRd_a][mo_range[0]:mo_range[1]]
        energies = molecule.electronic_structure.energies[mo_range[0]:mo_range[1]]
    else:
        ird_a_max = [np.argmax(ird_a_orb) for ird_a_orb in wfnsym.mo_IRd_a]
        energies = molecule.electronic_structure.energies

    ax1 = plt.axes()
    ax1.axes.get_xaxis().set_visible(False)  # Hide x axis
    # ax1.axes.get_yaxis().set_visible(True)

    degeneracy = [[energies[0]]]
    for energy in energies[1:]:
        if abs(energy - degeneracy[-1][-1]) < 1e-3:
            degeneracy[-1].append(energy)
        else:
            degeneracy.append([energy])

    max_value = 5e-3
    x_center = []
    for ix in degeneracy:
        if len(ix) == 1:
            x_center.append([0])
        else:
            x_center.append(np.linspace(-max_value, max_value, len(ix)))
    x_center = [y for x in x_center for y in x]

    plt.scatter(x_center, energies, s=500, marker="_", linewidth=3)
    for i in range(len(energies)):
        plt.text(-max_value * 2, energies[i], labels[ird_a_max[i]])

    plt.show()


def swap_vectors(v1, v2, position):
    vector1 = v1.copy()
    vector2 = v2.copy()
    for i in range(len(v1)):
        if i >= position:
            vector1[i] = v2[i]
            vector2[i] = v1[i]
    return vector1, vector2


def plot_symmetry_energy_evolution(molecules, wfnsym, mo_range=None):
    import matplotlib.pyplot as plt

    energies = []
    ird_a_max = []
    for idm, molecule in enumerate(molecules):
        labels = wfnsym[idm].IRLab
        if mo_range is not None:
            ird_a_max.append(np.array([np.argmax(ird_a_orb) for ird_a_orb in wfnsym[idm].mo_IRd_a]
                                      [mo_range[0]:mo_range[1]]))
            energies.append(molecule.electronic_structure.energies[mo_range[0]:mo_range[1]])
        else:
            ird_a_max.append(np.array([np.argmax(ird_a_orb) for ird_a_orb in wfnsym[idm].mo_IRd_a]))
            energies.append(molecule.electronic_structure.energies)

    energies_x_orbital = np.array(energies).T
    ird_a_x_orbital = np.array(ird_a_max).T

    for i in range(len(ird_a_x_orbital)):
        for j in range(len(ird_a_x_orbital[i])):
            if j == 0:
                old_ird = ird_a_x_orbital[i][0]
            else:
                if old_ird != ird_a_x_orbital[i][j]:
                    for k in range(len(ird_a_x_orbital) - i):
                        if old_ird == ird_a_x_orbital[k + i][j]:
                            ird_a_x_orbital[i], ird_a_x_orbital[k + i] = swap_vectors(ird_a_x_orbital[i],
                                                                                      ird_a_x_orbital[k + i], j)
                            energies_x_orbital[i], energies_x_orbital[k + i] = swap_vectors(energies_x_orbital[i],
                                                                                            energies_x_orbital[k + i],
                                                                                            j)
                            break
            old_ird = ird_a_x_orbital[i][j]

    for ide, energy in enumerate(energies_x_orbital):
        x = np.arange(len(energy))
        plt.plot(x, energy, marker='_')
        for i in range(len(energy)):
            plt.text(x[i], energy[i] + abs(energy[i])*0.001, labels[ird_a_x_orbital[ide][i]])

    plt.show()
