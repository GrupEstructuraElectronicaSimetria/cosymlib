import sys
from cosymlib import tools


def get_symgroup_data_txt(label, geometries, symgroup_results):
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


def get_symgroup_measure_txt(label, geometries, symgroup_results):

    sym_txt = 'Evaluating symmetry operation : {}\n \n'.format(label)
    for idx, geometry in enumerate(geometries):
        csm = symgroup_results[idx].csm
        max_name = len(max(geometry.get_name(), key=len))
        sym_txt += '{} '.format(geometry.get_name())
        if max_name < 9:
            n = 18 - len(geometry.get_name())
        else:
            n = 9 + max_name - len(geometry.get_name())
        sym_txt += '{:{width}.{prec}f}\n'.format(csm, width=n, prec=3)

    return sym_txt


def get_operated_matrices_txt(group, molecule, parsed_data):

    geometry = molecule.geometry
    sym_lables = parsed_data.SymLab
    sym_matrices = parsed_data.SymMat

    sep_line = '--------------------------------------------\n'
    txt_sym = 'MEASURES OF THE SYMMETRY GROUP:   {}\n'.format(group)
    txt_sym += 'Basis: {}\n'.format(list(molecule.electronic_structure.basis.keys())[0])
    txt_sym += sep_line
    txt_sym += ' Atomic Coordinates (Angstroms)\n'
    txt_sym += sep_line
    for idn, array in enumerate(geometry.get_positions()):
        txt_sym += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                              array[0], array[1], array[2])
    txt_sym += sep_line
    for i, group in enumerate(sym_lables):
        txt_sym += '\n'
        txt_sym += '@@@ Operation {0}: {1}'.format(i + 1, group)
        txt_sym += '\nSymmetry Transformation matrix\n'
        for array in sym_matrices[i]:
            txt_sym += ' {:11.6f} {:11.6f} {:11.6f}\n'.format(array[0], array[1], array[2])
        txt_sym += '\n'
        txt_sym += 'Symmetry Transformed Atomic Coordinates (Angstroms)\n'

        for idn, array in enumerate(geometry.get_positions()):
            array2 = np.dot(array, sym_matrices[i].T)
            txt_sym += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(geometry.get_symbols()[idn],
                                                                  array2[0], array2[1], array2[2])

    return txt_sym


def get_overlap_analysis_txt(parsed_data):

    ir_labels = parsed_data.SymLab
    ideal_gt = parsed_data.ideal_gt
    mo_soevs_a = parsed_data.mo_SOEVs_a
    mo_soevs_b = parsed_data.mo_SOEVs_b
    wf_soevs_a = parsed_data.wf_SOEVs_a
    wf_soevs_b = parsed_data.wf_SOEVs_b
    wf_soevs = parsed_data.wf_SOEVs
    grim_coef = parsed_data.grim_coef
    csm_coef = parsed_data.csm_coef

    sep_line = '   -------------------------------------------------------------------------------' \
               '------------------------------------------------------------------------\n'

    txt_sym = '\nIdeal Group Table\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    for i, line in enumerate(ideal_gt):
        txt_sym += '{:4}'.format(ir_labels[i]) + '  '.join(['{:7.3f}'.format(s) for s in line])
        txt_sym += '\n'
    txt_sym += sep_line

    txt_sym += '\nAlpha MOs: Symmetry Overlap Expectation Values\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line

    for i, line in enumerate(mo_soevs_a):
        txt_sym += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
        txt_sym += '\n'

    txt_sym += '\nBeta MOs: Symmetry Overlap Expectation Values\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    for i, line in enumerate(mo_soevs_b):
        txt_sym += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
        txt_sym += '\n'

    txt_sym += '\nWaveFunction: Symmetry Overlap Expectation Values\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    txt_sym += 'a-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_soevs_a])
    txt_sym += '\n'
    txt_sym += 'b-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_soevs_b])
    txt_sym += '\n'
    txt_sym += 'WFN ' + '  '.join(['{:7.3f}'.format(s) for s in wf_soevs])
    txt_sym += '\n'

    txt_sym += '\nWaveFunction: CSM-like values\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line

    txt_sym += 'Grim' + '  '.join(['{:7.3f}'.format(s) for s in grim_coef])
    txt_sym += '\n'
    txt_sym += 'CSM ' + '  '.join(['{:7.3f}'.format(s) for s in csm_coef])
    txt_sym += '\n'

    return txt_sym


def get_ir_analysis_txt(parsed_data):

    ir_labels = parsed_data.IRLab
    mo_ir_alpha = parsed_data.mo_IRd_a
    mo_ir_beta = parsed_data.mo_IRd_b
    wf_ir_a = parsed_data.wf_IRd_a
    wf_ir_b = parsed_data.wf_IRd_b
    wf_ir = parsed_data.wf_IRd

    sep_line = '   ---------------------------------------------\n'

    txt_sym = '\nAlpha MOs: Irred. Rep. Decomposition\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    for i, line in enumerate(mo_ir_alpha):
        txt_sym += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
        txt_sym += '\n'

    txt_sym += '\nBeta MOs: Irred. Rep. Decomposition\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    for i, line in enumerate(mo_ir_beta):
        txt_sym += '{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line])
        txt_sym += '\n'

    txt_sym += '\nWaveFunction: Irred. Rep. Decomposition\n'
    txt_sym += sep_line
    txt_sym += '     ' + '  '.join(['{:^7}'.format(s) for s in ir_labels])
    txt_sym += '\n'
    txt_sym += sep_line
    txt_sym += 'a-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_ir_a])
    txt_sym += '\n'
    txt_sym += 'b-wf' + '  '.join(['{:7.3f}'.format(s) for s in wf_ir_b])
    txt_sym += '\n'
    txt_sym += 'WFN ' + '  '.join(['{:7.3f}'.format(s) for s in wf_ir])
    txt_sym += '\n'

    return txt_sym
