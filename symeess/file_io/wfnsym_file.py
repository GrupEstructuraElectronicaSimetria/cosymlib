import numpy as np


def build_symmetry_operated_matrices(group, molecule, parsed_data):

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


def build_symmetry_overlap_analysis(parsed_data):

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


def build_symmetry_ireducible_representation_analysis(parsed_data):

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
