import sys
from cosymlib import tools


def build_symgroup_data(label, geometries, symgroup_results):
    # if output_name is not None:
    #     output = open(output_name + '.zout', 'w')
    #     output2 = open(output_name + '.ztab', 'w')
    # else:
    #     output = sys.stdout
    #     output2 = sys.stdout

    sep_line = '..................................................\n'

    sym_txt = 'Evaluating symmetry operation : {}\n'.format(label)

    for idx, geometry in enumerate(geometries):
        sym_txt += '{}\n'.format(geometry.get_name())
        sym_txt += '\n'
        sym_txt += 'Centered Structure\n'
        sym_txt += sep_line
        center_mass = tools.center_mass(geometry.get_symbols(), geometry.get_positions())
        for idn, array in enumerate(geometry.get_positions()):
            array = array - center_mass
            sym_txt += '{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(geometry.get_symbols()[idn],
                                                                  array[0], array[1], array[2])
        sym_txt += sep_line

        sym_txt += 'Optimal permutation\n'
        for idn, permutation in enumerate(symgroup_results[idx].optimum_permutation):
            sym_txt += '{:2} {:2}\n'.format(idn + 1, permutation)
        sym_txt += '\n'

        sym_txt += 'Inverted structure\n'
        for idn, axis in enumerate(symgroup_results[idx].nearest_structure):
            sym_txt += '{:2} {:12.8f} {:12.8f} {:12.8f}\n'.format(geometry.get_symbols()[idn],
                                                                    axis[0], axis[1], axis[2])
        sym_txt += '\n'

        sym_txt += 'Reference axis\n'
        for array in symgroup_results[idx].reference_axis:
            sym_txt += '{:12.8f} {:12.8f} {:12.8f}\n'.format(array[0], array[1], array[2])
        sym_txt += '\n'

        sym_txt += 'Symmetry measure {:.5f}\n'.format(symgroup_results[idx].csm)
        sym_txt += sep_line

    return sym_txt


def build_symgroup_measure(label, geometries, symgroup_results):

    sym_txt = 'Evaluating symmetry operation : {}\n \n'.format(label)
    for idx, geometry in enumerate(geometries):
        csm = symgroup_results[idx].csm
        max_name = len(max(geometry.get_name(), key=len))
        sym_txt += '{}'.format(geometry.get_name())
        if max_name < 9:
            n = 18 - len(geometry.get_name())
        else:
            n = 9 + max_name - len(geometry.get_name())
        sym_txt += '{:{width}.{prec}f}\n'.format(csm, width=n, prec=3)

    return sym_txt