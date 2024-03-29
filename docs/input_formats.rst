.. highlight:: rst
.. _input_formats:


Input file formats
==================

The computation of continuous shape or symmetry measures for molecules requires the input of
information on the geometry and the electronic structure of molecules which may be provided
in different formats. While the information on the molecular geometry is mandatory, providing
data on the electronic structure is optional and if this is not present, **cosymlib** may be
used to generate an approximate electronic structure (extended Hückel molecular orbitals or
a promolecular electronic density) using the molecular structure.

At this moment **cosymlib** is able to load structural data from files written in the
following formats:

1. Cartesian coordinate files (.xyz file)
2. Protein Data Bank Format (.pdb file)
3. Atomic coordinates generated by CSD's ConQuest (.cor file)

Information about the electronic structure can be supplied as a Gaussian formatted checkpoint
file (.fchk file) which can be generated by different quantum chemistry programs besides Gaussian.
In particular, the **Huckelpy** program developed within the **cosymlib** project can be used
to provide directly an electronic electronic structure for **cosymlib** or to generate a .fchk file
that can be used afterwards as an input for the scripts included in the **cosymlib**.

More detailed information on the file formats above is given in the sections below and the
corresponding links.

xyz files
^^^^^^^^^
The xyz format (`<http://en.wikipedia.org/wiki/XYZ_file_format>`_) is possybly the simplest
format used in computational chemistry to indicate a molecular geometry. To do so  you need
just to specify the number of atoms, a descriptive label and the cartesian coordinates for each
atom in the molecule, according to the following block structure:
::

 <number of atoms>
 <label>
 <element> <X> <Y> <Z>
  ...
 <element> <X> <Y> <Z>

The  <element> label for each atom is usually the corresponding atomic symbol, but this is only
strictly necessary for calculations involving the electronic structure. For geometric
CShMs or CSMs you may use arbitrary labels suc as A, B, C, ... for each vertex in the structure.

The units are generally given in Å, but this is not strictly necessary if you
are just interested in geometric shape and symmetry measures. Calculations involving the
electronic structure assume, however, that your coordinates are given in Å.

Example .xyz file:
::

 7
 CoCl6
 Co    0.0   0.0  0.0
 Cl   -2.4   0.0  0.0
 Cl    2.4   0.0  0.0
 Cl    0.0   2.4  0.0
 Cl    0.0  -2.4  0.0
 Cl    0.0   0.0  2.4
 Cl    0.0   0.0 -2.4

An .xyz file may contain one or several structures (all with the same number of atoms) in which case
you only need to include a different block for each structure, without leaving any blank line between
consecutive structures:
::

 n
 label_s1
 atom_s1_1 x_1 y_1 z_1
 ...
 atom_s1_n x_n y_n z_n
 n
 label_s2
 atom_s2_1 x_1 y_1 z_1
 ...
 atom_s2_n x_n y_n z_n

 ...


cor files
^^^^^^^^^

A .cor file contains atomic coordinates for a molecule or molecular
fragment generated by CSD's ConQuest (`<https://www.ccdc.cam.ac.uk/solutions/csd-core/components/conquest/>`_).
To generate the .cor file a search must be carried out in ConQuest in which only the atoms
in the fragment to be analyzed are defined. A coordinates file is then generated with the results
of the search ("Export entries as...") using the "Orthogonal" and "Hit Fragment Only" options.

Example of a CSD's .cor file:
::

 DOSDAQ **FRAG** 1
 Pd1  16.04450   0.00000   0.00000   1555011
 N1H  17.11627  -1.46619  -0.92480   9655010
 N1B  17.11627   1.46619   0.92483   3655010
 N1J  14.97273  -1.46619  -0.92480   1155501
 N1   14.97273   1.46619   0.92483   1555010
 FUBWUU **FRAG** 1
 Pt1   6.35325   1.50775   0.00000   1555001
 N1E*  4.67365   2.63675   0.26075   6665002
 N2    7.85296   1.80990   1.32660   1555016
 N2E*  4.85354   1.20560  -1.32660   6665016
 N1    8.03285   0.37875  -0.26075   1555002
 LEGFOS **FRAG** 1
 Cr1   5.47182   6.64300   6.21762   1555001
 N2    3.90370   5.54292   7.02349   1555003
 N1B   7.04319   5.53495   5.40674   3555002
 N2B   3.90370   7.74308   7.02349   3555003
 N1    7.04319   7.75105   5.40674   1555002


pdb files
^^^^^^^^^
The Protein Data Bank (pdb) file format is a textual file format describing
the three-dimensional structures of molecules held in the Protein Data Bank.
The pdb format accordingly provides for description and annotation of protein
and nucleic acid structures including atomic coordinates, secondary structure
assignments, as well as atomic connectivity, although the only section relevant for
**cosymlib** is the block containing the atomic coordinates.


