from cosym.molecule import Molecule, Geometry
from cosym import file_io
from cosym.file_io import shape2file
from cosym.utils import get_shape_map, molecular_orbital_diagram, symmetry_energy_evolution
import sys


class Cosym:
    """
    Main class of cosym program that can perform all the jobs
    """

    def __init__(self, structures):

        self._molecules = []
        if isinstance(structures, list):
            for structure in structures:
                if isinstance(structure, Molecule):
                    self._molecules.append(structure)
                elif isinstance(structure, Geometry):
                    self._molecules.append(Molecule(structure))
                else:
                    raise AttributeError('Molecule object not found')
        else:
            if isinstance(structures, Molecule):
                self._molecules.append(structures)
            elif isinstance(structures, Geometry):
                self._molecules.append(Molecule(structures))
            else:
                raise AttributeError('Molecule object not found')

    def write_shape_measure_2file(self, shape_reference, central_atom=0, output_name=None):
        """
        Method that prints to file shape's measure
        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's measure in the output_name.tab file
        """

        molecules_name = [molecule.get_name() for molecule in self._molecules]
        if type(shape_reference[0]) is Geometry:
            shape_results_measures = [self.get_shape_measure(reference.get_positions(), 'measure', central_atom)
                                      for reference in shape_reference]
            shape2file.write_shape_measure_data(shape_results_measures, molecules_name,
                                                [reference.get_name() for reference in shape_reference],
                                                output_name=output_name)
        else:
            shape_results_measures = [self.get_shape_measure(reference, 'measure', central_atom)
                                      for reference in shape_reference]
            shape2file.write_shape_measure_data(shape_results_measures, molecules_name,
                                                shape_reference, output_name=output_name)

    def write_shape_structure_2file(self, shape_reference, central_atom=0, output_name=None):
        """
        Method that prints to file shape's structure
        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """
        initial_geometries = [molecule.geometry for molecule in self._molecules]
        molecules_name = [molecule.get_name() for molecule in self._molecules]

        if type(shape_reference[0]) is not str:
            shape_results_structures = [self.get_shape_measure(reference.get_positions(), 'structure', central_atom)
                                        for reference in shape_reference]
            shape_results_measures = [self.get_shape_measure(reference.get_positions(), 'measure', central_atom)
                                      for reference in shape_reference]

            shape2file.write_shape_measure_data(shape_results_measures, molecules_name,
                                                [reference.get_name() for reference in shape_reference],
                                                output_name=output_name)
            geometries = []
            for idl, reference in enumerate(shape_reference):
                for idm, molecule in enumerate(self._molecules):
                    geometries.append(Geometry(symbols=molecule.geometry.get_symbols(),
                                               positions=shape_results_structures[idl][idm],
                                               name=reference.get_name()))
            file_io.write_file_xyz(geometries, output_name=output_name)
        else:
            shape_results_structures = [self.get_shape_measure(reference, 'structure', central_atom)
                                        for reference in shape_reference]
            shape_results_measures = [self.get_shape_measure(reference, 'measure', central_atom)
                                      for reference in shape_reference]

            shape2file.write_shape_measure_data(shape_results_measures, molecules_name,
                                                shape_reference, output_name=output_name)
            geometries = []
            for idl, reference in enumerate(shape_reference):
                for idm, molecule in enumerate(self._molecules):
                    geometries.append(Geometry(symbols=molecule.geometry.get_symbols(),
                                               positions=shape_results_structures[idl][idm],
                                               name=reference))
            file_io.write_file_xyz(geometries, output_name=output_name)

    def write_path_parameters_2file(self, shape_label1, shape_label2, central_atom=0,
                                    maxdev=15, mindev=0, maxgco=101, mingco=0, output_name=None):

        csm, devpath, GenCoord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom,
                                                            maxdev=maxdev, mindev=mindev, maxgco=maxgco, mingco=mingco)
        names_order = [molecule.get_name() for molecule in self._molecules]
        shape2file.write_minimal_distortion_path_analysis(shape_label1, shape_label2, csm, devpath, GenCoord,
                                                          maxdev, mindev, mingco, maxgco, names_order,
                                                          output_name=output_name)

    def write_symgroup_measure_all_info(self, group, multi=1, central_atom=0, output_name=None):
        if output_name is not None:
            output = open(output_name + '.zout', 'w')
        else:
            output = sys.stdout

        results = self.get_symgroup_measure(group=group, multi=multi, central_atom=central_atom)
        txt = file_io.header()
        txt += file_io.symgroup_file.build_symgroup_data(group, [molecule.geometry for molecule in self._molecules],
                                                         results)
        output.write(txt)

    def write_symgroup_measure(self, group, multi=1, central_atom=0, output_name=None):
        if output_name is not None:
            output = open(output_name + '.ztab', 'w')
        else:
            output = sys.stdout

        results = self.get_symgroup_measure(group=group, multi=multi, central_atom=central_atom)
        txt = file_io.header()
        txt += file_io.symgroup_file.build_symgroup_measure(group, [molecule.geometry for molecule in self._molecules],
                                                            results)
        output.write(txt)

    def write_wnfsym_measure_2file(self, group, vector_axis1=None, vector_axis2=None, center=None, output_name=None,
                                   n_molecule=0):
        if output_name is not None:
            output = open(output_name + '.wout', 'w')
        else:
            output = sys.stdout

        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        txt = file_io.header()
        txt += file_io.wfnsym_file.build_symmetry_operated_matrices(group, self._molecules[n_molecule],
                                                                    wfnsym_results[n_molecule])
        txt += file_io.wfnsym_file.build_symmetry_overlap_analysis(wfnsym_results[n_molecule])
        txt += file_io.wfnsym_file.build_symmetry_ireducible_representation_analysis(wfnsym_results[n_molecule])
        output.write(txt)

    def write_wnfsym_sym_matrices_2file(self, group, vector_axis1=None, vector_axis2=None, center=None, output_name=None,
                                        n_molecule=0):
        if output_name is not None:
            output = open(output_name + '.wout', 'w')
        else:
            output = sys.stdout

        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        txt = file_io.header()
        txt += file_io.wfnsym_file.build_symmetry_operated_matrices(group, self._molecules[0], wfnsym_results[n_molecule])
        output.write(txt)

    def write_wnfsym_sym_ovelap_2file(self, group, vector_axis1=None, vector_axis2=None, center=None, output_name=None,
                                      n_molecule=0):
        if output_name is not None:
            output = open(output_name + '.wout', 'w')
        else:
            output = sys.stdout

        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        txt = file_io.header()
        txt += file_io.wfnsym_file.build_symmetry_overlap_analysis(wfnsym_results[n_molecule])
        output.write(txt)

    def write_wnfsym_ireducible_repr_2file(self, group, vector_axis1=None, vector_axis2=None, center=None,
                                           output_name=None, n_molecule=0):
        if output_name is not None:
            output = open(output_name + '.wout', 'w')
        else:
            output = sys.stdout

        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        txt = file_io.header()
        txt += file_io.wfnsym_file.build_symmetry_ireducible_representation_analysis(wfnsym_results[n_molecule])
        output.write(txt)

    def write_mo_diagram(self, group, vector_axis1=None, vector_axis2=None, center=None,n_molecule=0):
        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        molecular_orbital_diagram(self._molecules[0], wfnsym_results[n_molecule])

    def write_sym_energy_evolution(self, group, vector_axis1=None, vector_axis2=None, center=None):
        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        symmetry_energy_evolution(self._molecules, wfnsym_results)

    def get_shape_measure(self, label, kind, central_atom=0):
        get_measure = 'get_shape_' + kind
        return [getattr(molecule.geometry, get_measure)(label, central_atom=central_atom)
                for molecule in self._molecules]

    def get_molecule_path_deviation(self, shape_label1, shape_label2, central_atom=0):
        return [molecule.geometry.get_path_deviation(shape_label1, shape_label2, central_atom) for molecule
                          in self._molecules]

    def get_molecule_generalized_coord(self, shape_label1, shape_label2, central_atom=0):
        return [molecule.geometry.get_generalized_coordinate(shape_label1, shape_label2, central_atom)
                for molecule in self._molecules]

    def get_path_parameters(self, shape_label1, shape_label2, central_atom=0, maxdev=15, mindev=0,
                            maxgco=101, mingco=0):

        def filter_results(results, criteria, max, min):
            filter = [True if x <= max and x >= min else False for x in criteria]
            results = [i for indx, i in enumerate(results) if filter[indx] == True]
            return results

        csm = {shape_label1: self.get_shape_measure(shape_label1, 'measure', central_atom),
               shape_label2: self.get_shape_measure(shape_label2, 'measure', central_atom)}
        devpath = self.get_molecule_path_deviation(shape_label1, shape_label2, central_atom)
        generalized_coord = self.get_molecule_generalized_coord(shape_label1, shape_label2, central_atom)
        criteria = devpath
        devpath = filter_results(devpath, criteria, maxdev, mindev)
        generalized_coord = filter_results(generalized_coord, criteria, maxdev, mindev)
        criteria = generalized_coord
        devpath = filter_results(devpath, criteria, maxgco, mingco)
        generalized_coord = filter_results(generalized_coord, criteria, maxgco, mingco)
        return csm, devpath, generalized_coord

    def get_symgroup_measure(self, group, multi=1, central_atom=0):
        return [molecule.geometry.get_symmetry_measure(label=group, multi=multi, central_atom=central_atom) for
                   molecule in self._molecules]

    def get_wfnsym_measure(self, group, vector_axis1, vector_axis2, center):
        return [molecule.get_mo_symmetry(group,  vector_axis1=vector_axis1,  vector_axis2=vector_axis2,  center=center)
                for molecule in self._molecules]

    def write_minimum_distortion_path_shape_2file(self, shape_label1, shape_label2, central_atom=0,
                                                  num_points=20, show=False, output_name=None):
        path_parameters = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom)[0]
        path = get_shape_map(shape_label1, shape_label2, num_points)
        if show:
            import matplotlib.pyplot as plt
            plt.plot(path[0], path[1], 'k', linewidth=2.0)
            plt.scatter(path_parameters[shape_label1], path_parameters[shape_label2], linewidths=0.01)
            plt.xlabel(shape_label1)
            plt.ylabel(shape_label2)
            plt.show()
        else:
            shape2file.write_shape_map(shape_label1, shape_label2, path, output_name)

    def write_point_group(self, tol=0.01):
        return [molecule.get_pointgroup(tol) for molecule in self._molecules]