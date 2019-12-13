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
    def measure(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].csm

    # Function description
    def nearest_structure(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].nearest_structure

    def optimum_axis(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].optimum_axis

    def optimum_permutation(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].optimum_permutation

    def reference_axis(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].reference_axis

    def cms_multi(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].cms_multi

    def axis_multi(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash].axis_multi

    def _do_measure(self, label, central_atom, multi, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if not symbols:
            self._symbols = None
        self._results[hash] = Symgroupy(self._coordinates,
                                        group=label,
                                        labels=self._symbols,
                                        central_atom=central_atom,
                                        multi=multi)

    def results(self, label, central_atom=None, multi=1, symbols=True):

        hash = hashlib.md5('{}{}'.format(label, central_atom).encode()).hexdigest()
        if hash not in self._results:
            self._do_measure(label, central_atom=central_atom, multi=multi, symbols=symbols)

        return self._results[hash]