[![Build Status](https://app.travis-ci.com/GrupEstructuraElectronicaSimetria/cosymlib.svg?branch=master)](https://app.travis-ci.com/github/GrupEstructuraElectronicaSimetria/cosymlib)
[![Coverage Status](https://coveralls.io/repos/github/GrupEstructuraElectronicaSimetria/cosymlib/badge.svg?branch=master)](https://coveralls.io/github/GrupEstructuraElectronicaSimetria/cosymlib?branch=master)
[![PyPI version](https://badge.fury.io/py/cosymlib.svg)](https://badge.fury.io/py/cosymlib)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4925766.svg)](https://doi.org/10.5281/zenodo.4925766)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GrupEstructuraElectronicaSimetria/cosymlib/)

CoSymLib
========
Python library to compute continuous shape and symmetry measures of molecular structures  
Online manual: http://cosymlib.readthedocs.io


Main features
-------------
- Continuous Shape and Symmetry measures of molecular geometries
- Continuous Symmetry measures of electronic structure of molecules
- Support for common file types including XYZ, coor, pdb & fchk 
- Usage through python API or command line tools

Requirements
------------
 - Python 2.7/3.4 or higher
 - Matplotlib
 - Numpy
 - PyYaml
 - Symgroupy
 - Wfnsympy
 - Huckelpy
 - Blas & Lapack libraries
 - Fortran compiler
 
Examples using command line tools 
---------------------------------
### Shape measure

````bash
$ shape cf4.xyz -c 1 -l

----------------------------------------------------------------------
 COSYMLIB v0.10.5
 Electronic Structure & Symmetry Group
 Institut de Quimica Teorica i Computacional (IQTC)
 Universitat de Barcelona
----------------------------------------------------------------------

Available reference structures with 4 Vertices:

Label       Sym       Info

SP-4        D4h       Square
T-4         Td        Tetrahedron
SS-4        C2v       Seesaw
vTBPY-4     C3v       Axially vacant trigonal bipyramid



$ shape cf4.xyz -c 1 -m all

----------------------------------------------------------------------
 COSYMLIB v0.10.5
 Electronic Structure & Symmetry Group
 Institut de Quimica Teorica i Computacional (IQTC)
 Universitat de Barcelona
----------------------------------------------------------------------

Structure         SP-4       T-4        SS-4       vTBPY-4

cf4,             33.333,     0.000,     9.790,     3.573,


----------------------------------------------------------------------
                    End of calculation
----------------------------------------------------------------------
````

### Geometrical symmetry measure
````bash
$ gsym cf4.xyz -l

----------------------------------------------------------------------
 COSYMLIB v0.10.6
 Electronic Structure & Symmetry Group
 Institut de Quimica Teorica i Computacional (IQTC)
 Universitat de Barcelona
----------------------------------------------------------------------

Available symmetry groups

E      Identity Symmetry
Ci     Inversion Symmetry Group
Cs     Reflection Symmetry Group
Cn     Rotational Symmetry Group (n: rotation order)
Sn     Rotation-Reflection Symmetry Group (n: rotation-reflection order)



$ gsym cf4.xyz -m Ci

----------------------------------------------------------------------
 COSYMLIB v0.10.6
 Electronic Structure & Symmetry Group
 Institut de Quimica Teorica i Computacional (IQTC)
 Universitat de Barcelona
----------------------------------------------------------------------

Evaluating symmetry operation : Ci
 
cf4          16.667

----------------------------------------------------------------------
                    End of calculation
----------------------------------------------------------------------
````

### Wave function symmetry analysis
````bash 
$ mosym NF3.fchk -m C3 -mo

----------------------------------------------------------------------
 COSYMLIB v0.10.6
 Electronic Structure & Symmetry Group
 Institut de Quimica Teorica i Computacional (IQTC)
 Universitat de Barcelona
----------------------------------------------------------------------


Symmetry group : C3
Molecule : NF3
--------------------------------------
Atomic Coordinates (Angstroms)
--------------------------------------
N    -0.000022    0.000014    0.483415
F     1.225707   -0.219543   -0.125324
F    -0.802975   -0.951706   -0.125331
F    -0.422715    1.171239   -0.125335

WaveFunction: Irred. Rep. Decomposition
     ------------------
        A        E   
     ------------------
a-wf  1.000    0.000
b-wf  1.000    0.000
WFN   1.000    0.000

Symmetry axis orientation
center:   0.00000036   -0.00000023    0.00244913
axis  :   0.00001508   -0.00001685    1.00000000
axis2 :   0.00000000    1.00000000    0.00001685

Molecule : NF3
--------------------------------------
Atomic Coordinates (Angstroms)
--------------------------------------
N    -0.000022    0.000014    0.483415
F     1.225707   -0.219543   -0.125324
F    -0.802975   -0.951706   -0.125331
F    -0.422715    1.171239   -0.125335

Alpha MOs: Irred. Rep. Decomposition
     ------------------------------------
       occup    E(eV)     A        E   
     ------------------------------------
   1     1   -708.54    0.994    0.006
   2     1   -708.53    0.002    0.998
   3     1   -708.53    0.004    0.996
   4     1   -424.19    1.000    0.000
   5     1    -45.23    1.000    0.000
   6     1    -41.57    0.000    1.000
   7     1    -41.57    0.000    1.000
   8     1    -26.74    1.000    0.000
   9     1    -18.48    0.000    1.000
  10     1    -18.48    0.000    1.000
  11     1    -18.40    1.000    0.000
  12     1    -15.18    0.000    1.000
  13     1    -15.18    0.000    1.000
  14     1    -14.24    1.000    0.000
  15     1    -13.54    0.000    1.000
  16     1    -13.54    0.000    1.000
  17     1    -10.69    1.000    0.000
  18     0     12.21    0.000    1.000
  19     0     12.21    0.000    1.000
  20     0     12.52    1.000    0.000

Beta MOs: Irred. Rep. Decomposition
     ------------------------------------
       occup    E(eV)     A        E   
     ------------------------------------
   1     1   -708.54    0.994    0.006
   2     1   -708.53    0.002    0.998
   3     1   -708.53    0.004    0.996
   4     1   -424.19    1.000    0.000
   5     1    -45.23    1.000    0.000
   6     1    -41.57    0.000    1.000
   7     1    -41.57    0.000    1.000
   8     1    -26.74    1.000    0.000
   9     1    -18.48    0.000    1.000
  10     1    -18.48    0.000    1.000
  11     1    -18.40    1.000    0.000
  12     1    -15.18    0.000    1.000
  13     1    -15.18    0.000    1.000
  14     1    -14.24    1.000    0.000
  15     1    -13.54    0.000    1.000
  16     1    -13.54    0.000    1.000
  17     1    -10.69    1.000    0.000
  18     0     12.21    0.000    1.000
  19     0     12.21    0.000    1.000
  20     0     12.52    1.000    0.000

----------------------------------------------------------------------
                    End of calculation
----------------------------------------------------------------------
````

Example using Python API
------------------------

````python
from cosymlib import Molecule, Geometry


# Define geometry
geometry = Geometry(positions=[[ 0.0000,  0.0000,  0.0000],
                               [ 0.5288,  0.1610,  0.9359],
                               [ 0.2051,  0.8240, -0.6786],
                               [ 0.3345, -0.9314, -0.4496],
                               [-1.0685, -0.0537,  0.1921]],
                    symbols=['C', 'H', 'H', 'H', 'H'])

# Shape measure
shp_measure = geometry.get_shape_measure('T-4', central_atom=1)

# Geometrical symmetry measure
sym_geom_measure = geometry.get_symmetry_measure('C3', central_atom=1)

#Create molecule from geometry (generate electronic structure with Extended Hukel method)
molecule = Molecule(geometry)

# Wave function symmetry measure
wf_sym_measure = molecule.get_wf_symmetry('Td')

# Electronic density measure
dens_sym_measure = molecule.get_dens_symmetry('Td')
````


Contact info
------------
Abel Carreras  
abelcarreras83@gmail.com  
Donostia International Physics Center (DIPC)

Pere Alemany  
p.alemany@ub.edu  
Electronic Structure & Symmetry group  
Department of Materials Science and Physical Chemistry  
Institut de Química Teòrica i Computacional (IQTC-UB)  
University of Barcelona
