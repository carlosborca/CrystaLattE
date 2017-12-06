#!/usr/bin/env python

import Read_CIF
import psi4
from psi4.driver.wrapper_autofrag import auto_fragments

# ======================================================================
# Read a CIF file and generates a supercell.
print ("")
print ("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ")
print (" STEP 1. GENERATION OF A SUPERCELL FROM A CIF FILE.")
print ("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ")
print ("")

ReadCifIn = 'Benzene-138K.cif'  # CIF input file.
ReadCifOut = 'Benzene-138K.xyz' # XYZ output file.
ReadCifA = '2'                  # Number of replicas on x.
ReadCifB = '2'                  # Number of replicas on y.
ReadCifC = '2'                  # Number of replicas on z.

args = ['', '-i', ReadCifIn, '-o', ReadCifOut, '-b', ReadCifA, ReadCifB, ReadCifC]

Read_CIF.main(args)
# ======================================================================

# ======================================================================
# Take the supercell .xyz file
# And generate .xyz files with all possible monomers.
#
# Check a possible solution for this at:
# http://www.psicode.org/psi4manual/master/psithonmol.html#advanced-python
print ("")
print ("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ")
print (" STEP 2. EXTRACTION OF MONOMERS FROM THE SUPERCELL.")
print ("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ")
print ("")

with open(ReadCifOut) as fxyz:
    sxyz = fxyz.readlines()
    #print(sxyz[2:])

SuperCell = psi4.geometry('\n'.join(sxyz[2:]))
SuperCell.update_geometry()
print (SuperCell.print_out())

CellFrags = auto_fragments(molecule = SuperCell)
print (CellFrags.natom())
print (CellFrags.nfragments())
# ======================================================================

# ======================================================================
# Loop through all monomers and generate dimers with all other monomers.
#
# Filter dimers that are too distant apart.
#
# Filter out and keep count of all non-unique dimers, using the nuclear
# repulsion energy criteria.
# ======================================================================

# ======================================================================
# Loop through all dimers and generate trimers with all other monomers.
#
# Filter trimers that are too distant apart.
#
# Filter out and keep count of all non-unique trimers, using ArbAlign.
# ======================================================================

# .
# .
# .

# ======================================================================
# Run plesantly parallel PSI4 computations on all the final list of 
# monomers, dimers, trimers, etc.
#
# Multiply the resulting energy of each one by the degeneracy factor.
#
# Sum results to get a lattice energy.
# ======================================================================
