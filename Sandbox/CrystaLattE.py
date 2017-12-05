#!/usr/bin/env python

# Read a CIF file and generates a supercell.
import Read_CIF
args = ['', '-i', 'Benzene-138K.cif', '-o',  'Benzene-138K.xyz', '-b', '5', '5', '5']
Read_CIF.main(args)

# Take the supercell .xyz file
# And generate .xyz files with all possible monomers.

# Loop through all monomers and generate dimers with all other monomers.
# Filter dimers that are too distant apart.
# Filter out and keep count of all non-unique dimers, using the nuclear repulsion energy criteria.

# Loop through all dimers and generate trimers with all other monomers.
# Filter trimers that are too distant apart.
# Filter out and keep count of all non-unique trimers, using ArbAlign RMSD criteria.

# .
# .
# .

# Run plesantly parallel PSI4 computations on all the final list of monomers, dimers, trimers,
# multiply the resulting energy of each one by the degeneracy factor
# Sum results to get a lattice energy.
