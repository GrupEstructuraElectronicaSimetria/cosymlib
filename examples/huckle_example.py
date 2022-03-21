from cosymlib import Cosymlib
from cosymlib.simulation import ExtendedHuckel
from cosymlib.molecule import Molecule
from cosymlib.molecule.geometry import Geometry


# Define geometry
geometry = Geometry(positions=[[ 0.0000,  0.0000,  0.0000],
                               [ 0.5288,  0.1610,  0.9359],
                               [ 0.2051,  0.8240, -0.6786],
                               [ 0.3345, -0.9314, -0.4496],
                               [-1.0685, -0.0537,  0.1921]],
                    symbols=['C', 'H', 'H', 'H', 'H'],
                    name='Methane')

# Build Molecule with ExtendedHuckel electronic structure
huckel_ee = ExtendedHuckel(geometry)
print('alpha electrons: ', huckel_ee.alpha_electrons)
print('beta electrons: ', huckel_ee.beta_electrons)
molecule = Molecule(geometry=geometry, electronic_structure=huckel_ee)

# Create Cosymlib instance and compute properties
mol = Cosymlib(molecule)
mol.print_esym_mo_irreducible_repr(group='td')
