from wfnsympy import WfnSympy
from cosymlib import tools


def _get_key(group, vector_axis1, vector_axis2, center):
    group_key = group.lower()
    vec1_key = ' '.join(['{:10.8f}'.format(n) for n in vector_axis1]) if vector_axis1 is not None else None
    vec2_key = ' '.join(['{:10.8f}'.format(n) for n in vector_axis2]) if vector_axis2 is not None else None
    center_key = ' '.join(['{:10.8f}'.format(n) for n in center]) if center is not None else None
    return group_key, vec1_key, vec2_key, center_key


class WfnSym:

    def __init__(self, molecule):
        self._molecule = molecule
        self._wfnsym_dict = self._molecule.electronic_structure.basis
        self._results = {}

    def symmetry_overlap_analysis(self, group, vector_axis1, vector_axis2, center):
        key = _get_key(group, vector_axis1, vector_axis2, center)
        if key not in self._results:
            self._do_measure(group, vector_axis1, vector_axis2, center)
        return [self._results[key].ideal_gt, self._results[key].SymLab, self._results[key].mo_SOEVs_a,
                self._results[key].mo_SOEVs_b, self._results[key].wf_SOEVs_a, self._results[key].wf_SOEVs_b,
                self._results[key].wf_SOEVs, self._results[key].grim_coef, self._results[key].csm_coef]

    def symmetry_irreducible_representation_analysis(self, group, vector_axis1, vector_axis2, center):
        key = _get_key(group, vector_axis1, vector_axis2, center)
        if key not in self._results:
            self._do_measure(group, vector_axis1, vector_axis2, center)
        return [self._results[key].IRLab, self._results[key].mo_IRd_a, self._results[key].mo_IRd_b,
                self._results[key].wf_IRd_a, self._results[key].wf_IRd_b, self._results[key].wf_IRd]

    def symmetry_matrix(self, group, vector_axis1, vector_axis2, center):
        key = _get_key(group, vector_axis1, vector_axis2, center)
        if key not in self._results:
            self._do_measure(group, vector_axis1, vector_axis2, center)
        return self._results[key].SymMat

    def _do_measure(self, group, vector_axis1, vector_axis2, center):
        key = _get_key(group, vector_axis1, vector_axis2, center)
        self._results[key] = WfnSympy(coordinates=self._molecule.geometry.get_positions(),
                                      symbols=self._molecule.geometry.get_symbols(),
                                      basis=self._wfnsym_dict,
                                      center=center, VAxis=vector_axis1, VAxis2=vector_axis2,
                                      alpha_mo_coeff=self._molecule.electronic_structure.coefficients_a,
                                      beta_mo_coeff=self._molecule.electronic_structure.coefficients_b,
                                      charge=self._molecule.electronic_structure.charge,
                                      multiplicity=self._molecule.electronic_structure.multiplicity,
                                      group=group.upper(),
                                      valence_only=self._molecule.electronic_structure.valence_only)

    def results(self, group, vector_axis1, vector_axis2, center):
        key = _get_key(group, vector_axis1, vector_axis2, center)
        if key not in self._results:
            self._do_measure(group, vector_axis1, vector_axis2, center)

        self._results[key].print_overlap_mo_alpha()
        return self._results[key]
