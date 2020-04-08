from cosymlib.symmetry import wfnsym, symgroup
from wfnsympy import WfnSympy
from symgroupy import Symgroupy
from cosymlib import tools
import numpy as np


def _get_key_symgroup(label, center, central_atom, connectivity, multi):
    group_key = label.lower()
    center_key = ' '.join(['{:10.8f}'.format(n) for n in center]) if center is not None else None
    connectivity_key = np.array2string(connectivity, precision=10) if connectivity is not None else None
    multi_key = int(multi)
    central_atom_key = int(central_atom)
    return group_key, center_key, central_atom_key, connectivity_key, multi_key


def _get_key_wfnsym(group, vector_axis1, vector_axis2, center):
    group_key = group.lower()
    vec1_key = ' '.join(['{:10.8f}'.format(n) for n in vector_axis1]) if vector_axis1 is not None else None
    vec2_key = ' '.join(['{:10.8f}'.format(n) for n in vector_axis2]) if vector_axis2 is not None else None
    center_key = ' '.join(['{:10.8f}'.format(n) for n in center]) if center is not None else None
    return group_key, vec1_key, vec2_key, center_key


class Symmetry:
    def __init__(self,
                 structure,
                 central_atom=None):

        # Allow geometry or molecule to be imported instead of crude Cartesian coordinates
        try:
            self._coordinates = structure.get_positions()
            self._symbols = structure.get_symbols()
        except AttributeError:
            try:
                self._coordinates = structure.geometry.get_positions()
                self._symbols = structure.geometry.get_symbols()
            except AttributeError:
                self._coordinates = structure
                self._symbols = None

        self._central_atom = central_atom
        self._results = {}

        try:
            self._electronic_structure = structure.electronic_structure
        except AttributeError:
            self._electronic_structure = None

    def get_symgroup_results(self, label, multi, connectivity=None, center=None, central_atom=None):

        if central_atom is None:
            central_atom = self._central_atom

        key = _get_key_symgroup(label, center, central_atom, connectivity, multi)
        if key not in self._results:
            self._results[key] = Symgroupy(self._coordinates,
                                           group=label,
                                           labels=self._symbols,
                                           central_atom=central_atom,
                                           multi=multi,
                                           center=center)
        return self._results[key]

    def get_wfnsym_results(self, group, vector_axis1, vector_axis2, center):

        if self._electronic_structure is None:
            print('Electronic structure not found')
            exit()

        key = _get_key_wfnsym(group, vector_axis1, vector_axis2, center)

        if key not in self._results:
            self._results[key] = WfnSympy(coordinates=self._coordinates,
                                          symbols=self._symbols,
                                          basis=self._electronic_structure.basis,
                                          center=center, VAxis=vector_axis1, VAxis2=vector_axis2,
                                          alpha_mo_coeff=self._electronic_structure.coefficients_a,
                                          beta_mo_coeff=self._electronic_structure.coefficients_b,
                                          charge=self._electronic_structure.charge,
                                          multiplicity=self._electronic_structure.multiplicity,
                                          group=group.upper(),
                                          valence_only=self._electronic_structure.valence_only)
        return self._results[key]

    ##########################################
    #       Structure symmetry methods       #
    ##########################################

    def measure(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).csm

    def nearest_structure(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).nearest_structure

    def optimum_axis(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).optimum_axis

    def optimum_permutation(self, label, multi=1, symbols=True):
        return self.get_symgroup_results(label, multi=multi).optimum_permutation

    def reference_axis(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).reference_axis

    def cms_multi(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).cms_multi

    def axis_multi(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).axis_multi

    ##########################################
    #       Electronic symmetry methods      #
    ##########################################

    def symmetry_overlap_analysis(self, group, vector_axis1, vector_axis2, center):
        results = self.get_wfnsym_results(group, vector_axis1, vector_axis2, center)

        return [results.ideal_gt, results.SymLab, results.mo_SOEVs_a,
                results.mo_SOEVs_b, results.wf_SOEVs_a, results.wf_SOEVs_b,
                results.wf_SOEVs, results.grim_coef, results.csm_coef]

    def symmetry_irreducible_representation_analysis(self, group, vector_axis1, vector_axis2, center):
        results = self.get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        return [results.IRLab, results.mo_IRd_a, results.mo_IRd_b,
                results.wf_IRd_a, results.wf_IRd_b, results.wf_IRd]

    def symmetry_matrix(self, group, vector_axis1, vector_axis2, center):
        results = self.get_wfnsym_results(group, vector_axis1, vector_axis2, center)
        return results.SymMat
