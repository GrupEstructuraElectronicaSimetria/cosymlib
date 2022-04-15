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

Python API example
------------------

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
