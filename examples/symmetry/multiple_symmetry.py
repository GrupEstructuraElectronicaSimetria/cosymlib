# This tutorial will show you how to handle a shape measure with lots of structures and later on their plot in what we
# call shape maps.
# Import part
from cosymlib.file_io import get_geometry_from_file_xyz, read_generic_structure_file
from cosymlib import Cosymlib

# This example will show the user how to perform a symmetry calculation on the ethane conformation path between
# the eclipsed and the staggered forms. The first thing to do is to read the xyz file that contains the molecules and
# creates one or multiple Geometry objects which will storage all the structural information

structures = get_geometry_from_file_xyz('../Symgroup/ethane.xyz', read_multiple=True)

# Now we are going to create a Cosymlib object with all these Geometry objects to simplify the printing process

structure_set = Cosymlib(structures)

# Finally, it is as simple as call one of the print functions inside Cosymlib class

structure_set.print_geometric_symmetry_measure('i')
