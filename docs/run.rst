.. highlight:: rst

=================
Executing program
=================

This script run symeess from terminal

.. module:: scripts/run

To run the program from terminal just use ::

   $ python run.py input_file


Shape
#####

The list of flags to run shape in symeess are the following ::

  * -m *Shape measure of input structure with reference polyhedra*
  * -s *Calculate the ideal structure for the input structure*
  * -t *Print the reference structure of the given label*

This flags requires the "-label reference_polyhedra".
Optional flags for shape are ::

  * -c int
  * -o output_name

where *int* is the position of the central atom in your molecule if exist.

Examples
********

If you want to make a simple measure of a tetrahedron molecule and compared to the ideal tetrahedron structure ::

  $ python run.py input_file -m -label T-4

Also you can do the same but for a tetrahedron molecule that have a central atom, i.e. in the first position ::

  $ python run.py input_file -m -label T-4 -c 1

To save it in a specific folder ::

  $ python run.py input_file -m -label T-4 -c 1 -o example_output_name

the output will be in the :file:`results` folder.
