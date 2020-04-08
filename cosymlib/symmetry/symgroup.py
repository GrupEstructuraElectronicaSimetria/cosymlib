from symgroupy import Symgroupy


def _get_key(label, center, central_atom, ignore_connectivity, multi):
    group_key = label.lower()
    center_key = ' '.join(['{:10.8f}'.format(n) for n in center]) if center is not None else None
    ignore_connectivity_key = bool(ignore_connectivity)
    multi_key = int(multi)
    central_atom_key = int(central_atom)
    return group_key, center_key, central_atom_key, ignore_connectivity_key, multi_key


class SymgroupXXX:

    def __init__(self, geometry, central_atom=None, connectivity=None):

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

        self._central_atom = central_atom
        self._connectivity = connectivity
        self._results = {}

    def get_symgroup_results(self, label, multi, ignore_connectivity=True, center=None, central_atom=None):

        if central_atom is None:
            central_atom = self._central_atom

        key = _get_key(label, center, central_atom, ignore_connectivity, multi)
        if key not in self._results:
            self._results[key] = Symgroupy(self._coordinates,
                                           group=label,
                                           labels=self._symbols,
                                           central_atom=central_atom,
                                           multi=multi,
                                           connectivity=self._connectivity,
                                           center=center)
        return self._results[key]

    # Function description
    def measure(self, label, multi=1):
        return self.get_symgroup_results(label, multi=multi).csm

    # Function description
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
