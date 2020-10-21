# This tutorial will show you how to handle a shape measure with lots of structures and later on their plot in what we
# call shape maps.
# Import part
from cosymlib.file_io import get_geometry_from_file_xyz, read_generic_structure_file
from cosymlib import Cosymlib

# In this case we got a scan of the dihedral angle of CrN6 molecule. The Oh and D3h shapes are going to be measure
# to analyse the path of the molecule along both shapes. In this case the structures will be read from an xyz file
# containing all the geometries along the path. This can be done in two ways, using the get_geometry_from_file_xyz
# function, which is specific for the xyz files or a more generic one, the read_generic_structure_file, that allows
# the user to read any kind of file supported by this program. In both cases we must change read_multiple to True,
# if not the function will only read the first structure in file

# structures = read_generic_structure_file('crn6.xyz', read_multiple=True)
structures = get_geometry_from_file_xyz('crn6.xyz', read_multiple=True)

# Creating shape maps is something that can be tedious to both new and veteran users. To avoid such task one can use the
# Cosymlib class, a class that allows the user to print all possible information that the program contains in a fancy
# format and also permits to plot shape maps more user friendly.

structure_set = Cosymlib(structures)

# Once a Cosymlib object is easy to print all the information in one map using the following function.

structure_set.print_minimum_distortion_path_shape('OC-6',
                                                  'TPR-6',
                                                  central_atom=1,
                                                  max_dev=30,
                                                  max_gco=105)
