from wfnsympy import WfnSympy
from symeess import tools
import hashlib


class Wfnsym:

    def __init__(self, molecule):

        self._molecule = molecule
        self._Ne_val = self._get_valence_electrons()
        self._wfnsym_dict = self._molecule.electronic_structure.basis
        self._results = {}

    def symmetry_overlap_analysis(self,
                                  label,
                                  vector_axis2,
                                  vector_axis1=[0., 0., 1.],
                                  center=[0., 0., 0.]):

        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, vector_axis1, vector_axis2, center)
        return [self._results[hash].ideal_gt, self._results[hash].SymLab, self._results[hash].mo_SOEVs_a,
                self._results[hash].mo_SOEVs_b, self._results[hash].mo_SOEVs, self._results[hash].wf_SOEVs_a,
                self._results[hash].wf_SOEVs_b, self._results[hash].wf_SOEVs, self._results[hash].grim_coef,
                self._results[hash].csm_coef]

    def symmetry_ireducible_representation_analysis(self, label,
                                                    vector_axis2,
                                                    vector_axis1=[0., 0., 1.],
                                                    center=[0., 0., 0.]):

        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, vector_axis1, vector_axis2, center)
        return [self._results[hash].IRLab, self._results[hash].mo_IRd_a, self._results[hash].mo_IRd_b,
                self._results[hash].wf_IRd_a, self._results[hash].wf_IRd_b, self._results[hash].wf_IRd]

    def symmetry_matrix(self, label,
                        vector_axis2,
                        vector_axis1=[0., 0., 1.],
                        center=[0., 0., 0.]):

        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, vector_axis1, vector_axis2, center)
        return self._results[hash].SymMat

    def _do_measure(self, label, vector_axis1, vector_axis2, center):

        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        self._results[hash] = WfnSympy(coordinates=self._molecule.geometry.get_positions(),
                                       symbols=self._molecule.geometry.get_symbols(),
                                       basis=self._wfnsym_dict,
                                       center=center, VAxis=vector_axis1, VAxis2=vector_axis2,
                                       alpha_mo_coeff=self._molecule.electronic_structure.coefficients_a,
                                       beta_mo_coeff=self._molecule.electronic_structure.coefficients_b,
                                       charge=self._molecule.electronic_structure.charge,
                                       multiplicity=self._molecule.electronic_structure.multiplicity,
                                       group=label.upper())

    def results(self, label, vector_axis2, vector_axis1=[0., 0., 1.], center=[0., 0., 0.]):

        hash = hashlib.md5('{}{}'.format(label, vector_axis1, vector_axis2, center).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, vector_axis1, vector_axis2, center)
        return self._results[hash]

    def _get_valence_electrons(self):
        n_valence = 0
        for symbol in self._molecule.geometry.get_symbols():
            n_valence += tools.element_valence_electron(symbol)
        return n_valence
