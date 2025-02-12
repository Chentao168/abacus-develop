from ase.build import fcc111
from ase.io import write

# Create the 2-layer Al(111) surface
atoms = fcc111('Al', size=(1, 1, 2))

# Add vacuum in the z-direction
atoms.center(vacuum=4.0, axis=2)

# Write the configuration to a file in POSCAR format
write('STRU', atoms, format='abacus')