.. highlight:: rst

How to use cosymlib
*******************

**Cosymlib** is a a python library for computing continuous symmetry & shape measures (CSMs & CShMs).
Besides using the APIs contained in **cosymlib** to build your own custom-made python programs we have
also written some general scripts to perform standard tasks such as calculating a continuous shape
measure for a given structure without need of writing a python script. All this general task scripts
are called using a similar syntax wich inludes the name of the script, the name of the file containing
the structural data and optional arguments specifying the tasks we want to perform::

   $ script filename -task1 -task2 ... -taskn

For instance, consider a ``struct.xyz`` file containing the following structural information
for a H\ :sub:`4`\  molecule in an approximately square geometry:

::

    4
    H4 Quadrangle
    H    1.1   0.9  0.0
    H   -1.0   1.1  0.0
    H   -0.9  -1.2  0.0
    H    1.1  -1.0  0.0

if we would like to compute the square shape measure S(SP-4) for this 4-vertex
polygon we simply can call the shape script indicating the name of the .xyz file
containing the coordinates and use the ``-mesure`` flag with the SP-4 label to indicate
that we want to use a square (SP-4 stands for square planar structure with 4 vertices)
as the reference shape::

   $ shape struct.xyz -measure SP-4

and the shape script wil call the APIs in cosymlib to read first your input file, generate a
molecule object, calculate the S(SP-4) continuous shape measure for it, and print
the result of the calculation on the screen.

If, for instance, we also want the coordinates for the square with the best overlap
with our problem structure, we just need to include the ``-structure`` option in our call::

   $ shape struct.xyz -measure SP-4 -structure

The general task scripts include also ``symmg`` and ``cchir`` for
calculating continuous symmetry and chirality measures for polyhedral structures, as well as
the ``wfsymm`` script for the continuous symmetry analysis wavefunctions and electron densities.

Besides these four basic scripts we have also included ``cosym``, a general script that allows to perform
any of the basic calculations above. We could, for instance, use directly ``cosym`` to calculate the previous
shape measure using the following command::

$ cosym struct.xyz -shape_measure SP-4 -structure

Notice that when using ``cosym`` some of the optional arguments in ``shape`` change to indicate which type
of calculation we would like to perform. For instance, ``-measure`` becomes ``-shape_measure``.
On the other hand, other arguments such as ``-structure``, which have also the same meaning when
calculating shape, symmetry or chirality measures, remain the same when used in combination with the
general ``cosym`` script.

Taking into account users of our previous programs, we have also written a stand-alone script
``shape_classic`` which is able to read an ``old_shape.dat`` input file containing both
the structural information and the necessary keywords to run a full CShM calculation as it was done
in our previous ``SHAPE`` program.

In the sections below you can find a detailed description of all stand-alone scripts as well as all APIs
included in the present distribution of **cosymlib**.

--------

General task scripts
--------------------
qqqqq

--------

cosym
^^^^^
qqqqq

--------

shape
^^^^^
qqqqq

--------

shape_classic
^^^^^^^^^^^^^
To run ``shape_classic`` you only need an ``old_shape.dat`` input file containing both
the structural information and the necessary keywords to run a full CShM calculation as in
the old ``SHAPE`` program::

   $ shape_classic old_shape.dat

The script will perform all tasks indicated in the input file, creating the necessary output
files, normally ``old_shape.out`` and ``old_shape.tab`` with the same information as when using
our previous ``SHAPE`` program.  Follow the link below  for a pdf version of the user guide for
SHAPE ver. 2.1 where you will find all information to perform a continuous shape analysis using
this option.

:download: `SHAPE ver. 2.1 User's guide.  <User_Manual_SHAPE_v2.1.pdf>`



--------

symmg
^^^^^
qqqqq

--------

cchir
^^^^^
qqqqq

--------

wfsymm
^^^^^^
qqqqq

--------

Specific task scripts
---------------------

qqqqq

--------

shape_map
^^^^^^^^^
qqqqq

--------

Using cosymlib's APIs
---------------------

qqqqq

--------



