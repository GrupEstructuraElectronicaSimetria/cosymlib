# This tutorial will show you how to handle a simple shape measure using a python script. It will not cover the
# reading and wrinting file which is explined in detail in the ----
# Import part:
from cosymlib import Geometry

# In this case we got a fac-[FeCl3Br3]3- slightly distorted from his octahedral shape and we want to know how distorted
# it is with continious shape measure. As we are analyzing shapes and not symmetry, we can avoid chemical part
# (elements and electronic structure).

structure = [[-0.000250,  -0.000142,  -0.002855],
             [ 2.028556,   0.001959,   1.468881],
             [-1.015070,  -1.758645,   1.467680],
             [-1.016235,   1.757020,   1.467903],
             [-2.158748,   0.000226,  -1.487093],
             [ 1.079740,   1.869427,  -1.486814],
             [ 1.082006,  -1.869845,  -1.484507]]
symbols = ['Fe', 'Cl', 'Cl', 'Cl', 'Br', 'Br', 'Br']

# The first thing we should do is to create a Geometry or Molecule object. As stated before, as we do not need the
# chemical part, the geometry should be our choice, as it only focus on the structural part.
# The Geometry class, only needs the position of the atoms in order to be created.

fecl3br3 = Geometry(structure, symbols=symbols)

# Now, different approaches can be performed, but as this tutorial try to be as simple as possible, the most striaght
# approach is to just call for the shape measure within the Geometry class. Shape only requires an ideal or a user's
# reference structure for the comparison with the case structure. However, as this molecule has a central atom (Fe3+),
# it is a must to include this information in the function call.

measure = fecl3br3.get_shape_measure('OC-6', central_atom=1)
print('Oh measure: {:.2f}'.format(measure))
