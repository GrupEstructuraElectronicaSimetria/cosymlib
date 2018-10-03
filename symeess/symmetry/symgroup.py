import numpy as np
import hashlib
from symgroupy import Symgroupy


class Symgroup:

    def __init__(self, geometry):

        # Allow geometry or molecule to be imported instead of crude Cartesian coordinates
        try:
            self._coordinates = geometry.get_positions()
            self._symbols = geometry.get_symbols()
        except AttributeError:
            try:
                self._coordinates = geometry.geometry.get_positions()
                self._symbols = geometry.geometry.get_symbols()
            except AttributeError:
                self._coordinates = geometry
                self._symbols = None

        self._coordinates = np.ascontiguousarray(self._coordinates)

        self._results = {}

    # Function description
    def get_measure(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].csm

    # Function description
    def get_nearest_structure(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_nearest_structure

    def get_optimum_axis(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_optimum_axis

    def get_optimum_permutation(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_optimum_permutation

    def get_reference_axis(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_reference_axis

    def get_cms_multi(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_cms_multi

    def get_axis_multi(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash].get_axis_multi

    def _do_measure(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        self._results[hash] = Symgroupy(self._coordinates,
                                        group=label,
                                        labels=self._symbols,
                                        central_atom=central_atom,
                                        multi=multi)

    def get_results(self, label, central_atom=None, multi=None):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi)

        return self._results[hash]