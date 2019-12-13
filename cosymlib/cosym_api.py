from cosymlib.molecule import Molecule, Geometry
from cosymlib import file_io
from cosymlib.file_io import shape2file
from cosymlib.utils import get_shape_map, molecular_orbital_diagram, symmetry_energy_evolution
import sys


class Cosymlib:
    """
    Main class of cosymlib program that can perform all the jobs
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

    def write_shape_measure_2file(self, shape_reference, central_atom=0, fix_permutation=False, output_name=None):
        """
        Method that prints to file shape's measure
        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's measure in the output_name.tab file
        """

        if output_name is not None:
            output = open(output_name + '_tab.csv', 'w')
        else:
            output = sys.stdout

        molecules_name = [molecule.get_name() for molecule in self._molecules]
        shape_results_measures = []
        references = []
        for reference in shape_reference:
            if type(reference) is Geometry:
                shape_results_measures.append(self.get_shape_measure(reference.get_positions(), 'measure', central_atom,
                                                                     fix_permutation))
                references.append(reference.get_name())
            else:
                shape_results_measures.append(self.get_shape_measure(reference, 'measure', central_atom,
                                                                     fix_permutation))
                references.append(reference)

        output.write(file_io.header())
        output.write(shape2file.write_shape_measure_data(shape_results_measures, molecules_name, references))

    def write_shape_structure_2file(self, shape_reference, central_atom=0, fix_permutation=False, output_name=None):
        """
        Method that prints to file shape's structure
        :param shape_reference: reference polyhedra label which user will compare with his polyhedra.
                            Reference labels can be found in [#f1]_
        :param central_atom: position of the central atom in molecule if exist
        :param output_name: custom name without extension
        :return: shape's structure in the output_name.out file
        """

        if output_name is not None:
            output = open(output_name + '_shp.xyz', 'w')
        else:
            output = sys.stdout

        shape_results_structures = []
        references = []
        for reference in shape_reference:
            if type(reference) is Geometry:
                shape_results_structures.append(self.get_shape_measure(reference.get_positions(),
                                                                       'structure', central_atom, fix_permutation))
                references.append(reference.get_name())
            else:
                shape_results_structures.append(self.get_shape_measure(reference, 'structure', central_atom,
                                                                       fix_permutation))
                references.append(reference)

        self.write_shape_measure_2file(shape_reference, central_atom, fix_permutation, output_name)
        geometries = []
        for idm, molecule in enumerate(self._molecules):
            geometries.append(molecule.geometry)
            for idl, reference in enumerate(references):
                geometries.append(Geometry(symbols=molecule.geometry.get_symbols(),
                                           positions=shape_results_structures[idl][idm],
                                           name=molecule.get_name() + '_' + reference))
        output.write(file_io.write_file_xyz(geometries))

    def write_path_parameters_2file(self, shape_label1, shape_label2, central_atom=0,
                                    maxdev=15, mindev=0, maxgco=101, mingco=0, output_name=None):

        if output_name is not None:
            output = open(output_name + '_tab.csv', 'w')
        else:
            output = sys.stdout

        output.write(file_io.header())
        csm, devpath, GenCoord = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom,
                                                          maxdev=maxdev, mindev=mindev, maxgco=maxgco, mingco=mingco)
        names_order = [molecule.get_name() for molecule in self._molecules]
        txt = shape2file.write_minimal_distortion_path_analysis(csm, devpath, GenCoord, maxdev, mindev,
                                                                mingco, maxgco, names_order)
        output.write(txt)

    def write_symgroup_measure_all_info(self, group, multi=1, central_atom=0, symbols=True, output_name=None):

        if output_name is not None:
            output = open(output_name + '.zout', 'w')
        else:
            output = sys.stdout

        results = self.get_symgroup_measure(group=group, multi=multi, symbols=symbols, central_atom=central_atom)
        txt = file_io.header()
        txt += file_io.symgroup_file.build_symgroup_data(group, [molecule.geometry for molecule in self._molecules],
                                                         results)
        output.write(txt)

    def write_symgroup_measure(self, group, multi=1, central_atom=0, symbols=True, output_name=None):

        if output_name is not None:
            output = open(output_name + '.ztab', 'w')
        else:
            output = sys.stdout

        results = self.get_symgroup_measure(group=group, multi=multi, symbols=symbols, central_atom=central_atom)
        txt = file_io.header()
        txt += file_io.symgroup_file.build_symgroup_measure(group, [molecule.geometry for molecule in self._molecules],
                                                            results)
        output.write(txt)

    def write_symgroup_structure(self, group, multi=1, central_atom=0, symbols=True, output_name=None):

        if output_name is not None:
            output = open(output_name + '_sym.xyz', 'w')
        else:
            output = sys.stdout

        self.write_symgroup_measure(group=group, multi=multi, symbols=symbols, central_atom=central_atom)
        results = self.get_symgroup_measure(group=group, multi=multi, symbols=symbols, central_atom=central_atom)
        geometries = []
        for idm, molecule in enumerate(self._molecules):
            geometries.append(Geometry(symbols=molecule.geometry.get_symbols(),
                                       positions=results[idm].nearest_structure,
                                       name=molecule.get_name() + '_' + group + ' with orientation ' +
                                            ' '.join('{:.8f}'.format(e) for e in results[idm].optimum_axis)))
        output.write(file_io.write_file_xyz(geometries))

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

    def write_wnfsym_sym_matrices_2file(self, group, vector_axis1=None, vector_axis2=None, center=None,
                                        output_name=None,
                                        n_molecule=0):
        if output_name is not None:
            output = open(output_name + '.wout', 'w')
        else:
            output = sys.stdout

        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        txt = file_io.header()
        txt += file_io.wfnsym_file.build_symmetry_operated_matrices(group, self._molecules[0],
                                                                    wfnsym_results[n_molecule])
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

    def write_mo_diagram(self, group, vector_axis1=None, vector_axis2=None, center=None, n_molecule=0):
        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        molecular_orbital_diagram(self._molecules[0], wfnsym_results[n_molecule])

    def write_sym_energy_evolution(self, group, vector_axis1=None, vector_axis2=None, center=None):
        wfnsym_results = self.get_wfnsym_measure(group, vector_axis1, vector_axis2, center)
        symmetry_energy_evolution(self._molecules, wfnsym_results)

    def get_shape_measure(self, label, kind, central_atom=0, fix_permutation=False):
        get_measure = 'get_shape_' + kind
        return [getattr(molecule.geometry, get_measure)(label, central_atom=central_atom,
                                                        fix_permutation=fix_permutation)
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

        if type(shape_label1) is Geometry:
            label1 = shape_label1.get_positions()
            label1_name = shape_label1.get_name()
        else:
            label1 = shape_label1
            label1_name = shape_label1
        if type(shape_label2) is Geometry:
            label2 = shape_label2.get_positions()
            label2_name = shape_label2.get_name()
        else:
            label2 = shape_label2
            label2_name = shape_label2

        csm = {label1_name: self.get_shape_measure(label1, 'measure', central_atom),
               label2_name: self.get_shape_measure(label2, 'measure', central_atom)}
        devpath = self.get_molecule_path_deviation(label1, label2, central_atom)
        generalized_coord = self.get_molecule_generalized_coord(label1, label2, central_atom)
        criteria = devpath
        devpath = filter_results(devpath, criteria, maxdev, mindev)
        generalized_coord = filter_results(generalized_coord, criteria, maxdev, mindev)
        criteria = generalized_coord
        devpath = filter_results(devpath, criteria, maxgco, mingco)
        generalized_coord = filter_results(generalized_coord, criteria, maxgco, mingco)
        return csm, devpath, generalized_coord

    def get_symgroup_measure(self, group, multi=1, central_atom=0, symbols=True):
        return [molecule.geometry.get_symmetry_measure(label=group, multi=multi, central_atom=central_atom,
                                                       symbols=symbols) for molecule in self._molecules]

    def get_wfnsym_measure(self, group, vector_axis1, vector_axis2, center):
        return [molecule.get_mo_symmetry(group, vector_axis1=vector_axis1, vector_axis2=vector_axis2, center=center)
                for molecule in self._molecules]

    def write_minimum_distortion_path_shape_2file(self, shape_label1, shape_label2, central_atom=0,
                                                  num_points=20, output_name=None):

        self.write_path_parameters_2file(shape_label1, shape_label2, central_atom=central_atom, output_name=output_name)
        if output_name is not None:
            output = open(output_name + '_pth.csv', 'w')
            output2 = open(output_name + '_pth.xyz', 'w')
        else:
            output = sys.stdout
            output2 = sys.stdout

        if type(shape_label1) is Geometry:
            label1 = shape_label1.get_positions()
            label1_name = shape_label1.get_name()
        else:
            label1 = shape_label1
            label1_name = shape_label1
        if type(shape_label2) is Geometry:
            label2 = shape_label2.get_positions()
            label2_name = shape_label2.get_name()
        else:
            label2 = shape_label2
            label2_name = shape_label2

        path_parameters = self.get_path_parameters(shape_label1, shape_label2, central_atom=central_atom)[0]
        path = get_shape_map(label1, label2, num_points)
        output.write(shape2file.write_shape_map(label1_name, label2_name, path))
        test_structures = []
        for ids, structure in enumerate(path[2]):
            test_structures.append(Geometry(symbols=['' for _ in range(len(structure))],
                                            positions=structure, name='map_structure{}'.format(ids)))
        output2.write(file_io.write_file_xyz(test_structures))
        if output_name is None:
            import matplotlib.pyplot as plt
            plt.plot(path[0], path[1], 'k', linewidth=2.0)
            plt.scatter(path_parameters[label1_name], path_parameters[label2_name], linewidths=0.01)
            plt.xlabel(label1_name)
            plt.ylabel(label2_name)
            plt.show()

    def write_point_group(self, tol=0.01):
        return [molecule.get_pointgroup(tol) for molecule in self._molecules]
