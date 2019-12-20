.. highlight:: rst

=================
Executing program
=================

This script run symeess from terminal

.. module:: scripts/cosym

To run the program from terminal just use ::

   $ symeess_script


Shape
#####

This flags requires the "-label reference_polyhedra" flag.
The list of flags to run shape in symeess are the following ::

  · -m *Shape measure of input structure with reference polyhedra*
  · -s *Calculate the ideal structure for the input structure*
  · -t *Print the reference structure of the given label*

Optional flags for shape are ::

  · -c int\ *where int is the postion of the central atom in your structure*
  · -o output_name *output name without extension*

Extra flags ::

  · -n int\ *generate the possible ideal structures with int vertices stored in the program*
  · -map *calculates the minimal distortion interconversion path between two given polyhedra*
  · -path *calculates the path deviation function as well as the generalized coordinate of the studied
  structures in a given path*

Examples
********

If you want to make a simple measure of a tetrahedron molecule and compared to the ideal tetrahedron structure ::

  $ symeess_script -m -label T-4 -input_file input_name

Also you can do the same but for a tetrahedron molecule that have a central atom, i.e. in the first position ::

  $ symeess_script -m -label T-4 -c 1 -input_file input_name

To save it with a specific name ::

  $ symeess_script -m -label T-4 -c 1 -output_name example_name -input_file input_name

To calculate the map between the tetrahedron and the square ::

  $ symeess_script -map -label T-4 SP-4

the output will be in the :file:`results` folder.

Symgroup
########

This flags requires the "-group reference_polyhedra" flag.
The list of flags to run shape in symeess are the following ::

  · -symm *Symgroup measure of input structure with reference polyhedra*

