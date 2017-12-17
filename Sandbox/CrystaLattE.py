#!/usr/bin/env python

import fnmatch
import math
import psi4
import re
import os
import sys
from psi4.driver.wrapper_autofrag import auto_fragments

# Imports of outsourced code.
sys.path.insert(0, "Read_CIF")
import Read_CIF

# ==================================================================
def print_error(msg):
    print("")
    print(" ERROR:  %s" % msg)
    print("")
# ==================================================================

# ==================================================================
def psifrg2xyz(p4frag, numfrags, numfatoms, frg_separator):
    """Takes the output of the automatic fragmentation procedure, the total number of fragments,
    the number of atoms in one fragment, and the string that separates the fragments.
    Returns .xyz files for each fragment."""

    # Read the output of the automatic fragmentation.
    with open(p4frag) as cellfrags:
        frags = cellfrags.readlines()

    # Generate .xyz files for each fragment.
    fcounter = 1         # Counter of processed fragments.
    lcounter = 0         # Counter of lines.
    HeaderLine = True    # Flag for first line of a fragment .xyz file.

    for line in frags:
        lcounter += 1

        if line.startswith(frg_separator):
            fcounter += 1
            HeaderLine = True

        else:
            ffnidx = "f" + str(fcounter).zfill(len(str(numfrags))) \
                     + ".xyz"
            with open(ffnidx, "a") as frgxyz: # WARNING: File exists?

                if HeaderLine == True:
                    frgxyz.write(str(numfatoms) + "\n" + "Fragment " # TODO: What if there are more than two molecules in the central unit cell.
                    + str(fcounter) + "\n")
                    HeaderLine = False

                if HeaderLine == False:
                    frgxyz.write(line)
    
    print("%i fragment .xyz files were created." % fcounter)

    return
# ==================================================================

# ==================================================================
def frags_filter(n_atm_frg): # WARNING: The number of atoms in a fragment probably should be passed as a list.
    """It takes the possible numbers of atoms in a fragment and discards all fragments files that contain a smaller number of atoms."""

    # TODO: What if there are more than two types of molecules in the
    #       unit cell? n_atom_frg should be passed as a list with
    #       possible numbers of atoms in a fragment.

    directory = os.getcwd()
    files = os.listdir(directory)
    incmolcounter = 0

    for file in files:

        if re.match('^f[0-9]+.xyz$', file): # Match filenames f???.xyz

            with open(file, 'r') as f:
                lines = f.readlines()

                if len(lines) != n_atm_frg + 2: # Match files with unexpected atoms.
                    print("Expected " + str(n_atm_frg)\
                          + " atoms. Found only " + str(len(lines) - 2)\
                          + ". Removing: " + file)

                    os.remove(os.path.join(directory,file))

                    incmolcounter += 1

    print("\nRemoved %s fragments.\n" % incmolcounter)

    return
# ==================================================================

# ==================================================================
def scellcntr(xyzpattern):
    """Compute the center of coordinates of .xyz files with a matching pattern in the filename.
    Returns a list with fragment names, x, y, and z coordinates in lists of floats and the coordinates of the center."""

    frg = [] # Fragments file names.
    x = [] # Values of x coordinates.
    y = [] # Values of y coordinates.
    z = [] # Values of z coordinates.

    directory = os.getcwd()
    files = os.listdir(directory)

    for file in files:

        if re.match(xyzpattern, file):
            
            with open(file, 'r') as f:
                lines = f.readlines()

                i = 0
                for i in range(len(lines)):

                    if i < 2:
                        continue

                    else:
                        info = lines[i].split()
                        frg.append(file)
                        x.append(float(info[1]))
                        y.append(float(info[2]))
                        z.append(float(info[3]))

    cntr_x = (max(x) - min(x))/2.0
    cntr_y = (max(y) - min(y))/2.0
    cntr_z = (max(z) - min(z))/2.0

    print("Center of the supercell located at:")
    print("x= %10.5f"   % (cntr_x))
    print("y= %10.5f"   % (cntr_y))
    print("z= %10.5f\n" % (cntr_z))

    return frg, x, y, z, cntr_x, cntr_y, cntr_z 
# ==================================================================

# ==================================================================
def posvecmag(xyzpattern):
    """Using the scellcntr() function it computes the magnitude of the position vector of the atoms read by the scellcntr() function.
    Returns the filename and the magnitude of the position vector."""

    frg, x, y, z, cntr_x, cntr_y, cntr_z = scellcntr(xyzpattern)

    r = []

    for j in range(len(x)):
        r.append(math.sqrt((x[j] - cntr_x)**2.0 + (y[j] - cntr_y)**2.0 + (z[j] - cntr_z)**2.0))

    print("Magnitudes of the position vectors calculated.\n")

    return frg, r
# ==================================================================

# ==================================================================
def proximity():
    """Using position vector magnitude, computed in posvecmag(). This function returns a list of fragments filenames
    with fragments containing atoms closest to the center of the supercell first."""

    frg, r = posvecmag("^f[0-9]+.xyz$")

    frgprox = []

    for pair in sorted(enumerate(r), key=lambda item:item[1]):
        if frg[pair[0]] in frgprox:
            continue
        else:
            frgprox.append(frg[pair[0]])

    print("Molecules organized based on the closest atom to the")
    print("center of the supercell.\n")

    return frgprox
# ==================================================================

# ==================================================================
def frgs2mols():
    """Takes fragents and a list of proximity to the center of the supercell, generated by the proximity() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell."""

    directory = os.getcwd()
    proxlist = proximity()
    molecule = 0

    for file in proxlist:
        molecule += 1
        molfname = "m" + str(molecule).zfill(len(proxlist[-1]) - 5) + ".xyz"

        with open(file, 'r') as ffile, open(molfname, 'w') as mfile:

            for l in ffile.readlines():

                if l.startswith("Fragment"):
                    newl = "Molecule " + str(molecule) + " (" + l[:-1] + ")\n"
                    mfile.write(newl)

                else:
                    mfile.write(l)

        os.remove(os.path.join(directory, file))

        print("Fragment " + str(file) + " has been rewritten to: " + str(molfname))

    print("")
# ==================================================================

# ==================================================================
def mols2mons():
    """Takes molecules and filters and a list of proximity to the center of the supercell, generated by the proximity() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell.
    It returns four global variables with the number of monomeers and the maximum possible number of dimers, trimers, tetramers."""

    directory = os.getcwd()
    molfiles = os.listdir(directory)
    monidx = 0
    
    global total_monomers
    total_monomers = 0

    for file in molfiles:
        
        if re.match('^m[0-9]+.xyz$', file):
            monidx += 1
            monfname = "1-" + str(file)[1:]

            with open(file, 'r') as moleculef, open(monfname, 'w') as monomerf:
                total_monomers += 1

                for l in moleculef.readlines():

                    if l.startswith("Molecule"):
                        newl = "Monomer " + str(monidx) + "\n"
                        monomerf.write(newl)

                    else:
                        monomerf.write(l)

            os.remove(os.path.join(directory, file))

            print("Molecule " + str(file) + " has been rewritten to: " + str(monfname))

    print("\nTotal number of mononmers: %i" % total_monomers)
    
    global total_dimers
    total_dimers = (total_monomers)*(total_monomers - 1)/2

    print("\nMaximum possible number of dimers: %i" % total_dimers)

    global total_trimers
    total_trimers = (total_monomers)*(total_monomers - 1)*(total_monomers - 2)/6

    print("\nMaximum possible number of trimers: %i" % total_trimers)

    global total_tetramers
    total_tetramers = (total_monomers)*(total_monomers - 1)*(total_monomers - 2)*(total_monomers - 3)/24

    print("\nMaximum possible number of tetramers: %i\n" % total_tetramers)

# ==================================================================

# ==================================================================
def rmin(nmer1, nmer2):
    """Takes two N-mer files and returns separation between closest atoms belonging to different N-mers."""

    with open(nmer1, 'r') as nmf1, open(nmer2, 'r') as nmf2:
        nm1l = nmf1.readlines()
        nm2l = nmf2.readlines()

        nmf1info = []
        nmf2info = []

        dist = []

        i = 0
        for i in range(len(nm1l)):

            if i < 2:
                continue

            else:
                nmf1atom = nm1l[i].split()
                nmf1info.append(nmf1atom)

        for j in range(len(nm2l)):

            if j < 2:
                continue

            else:
                nmf2atom = nm2l[j].split()
                nmf2info.append(nmf2atom)

        for atom1 in range(len(nmf1info)):

            for atom2 in range(len(nmf2info)):
                x12 = float(nmf1info[atom1][1]) - float(nmf2info[atom2][1])
                y12 = float(nmf1info[atom1][2]) - float(nmf2info[atom2][2])
                z12 = float(nmf1info[atom1][3]) - float(nmf2info[atom2][3])
                r12 = math.sqrt(x12**2.0 + y12**2.0 + z12**2.0)
                dist.append(r12)

    minsep = min(dist)

    return minsep
# ==================================================================

# ==================================================================
def nmermerger(nm_a, nm_b, nm_new):
    """Takes the .xyz filenames of two existing N-mers and of a new N-mer to be written.
    Produces the .xyz file of the new N-mer."""

    with open(nm_a, 'r') as nm_a_f, open(nm_b, 'r') as nm_b_f:
        
        nm_n_l = []

        i = 0
        for line_nm_a_f in nm_a_f.readlines():
            i += 1

            if i == 1:
                nm_n_l_line = str(int(line_nm_a_f) + int(nm_b_f.readline())) + "\n"
                nm_n_l.append(nm_n_l_line)

            elif line_nm_a_f.startswith("Monomer"):
                nm_n_l_new_line = "Monomers " + str(nm_a)[2:-4] + "+" + str(nm_b)[2:-4] + "\n"
                nm_n_l.append(nm_n_l_new_line)

            else:
                nm_n_l.append(line_nm_a_f)
        
        nm_n_l.append("\n")

        j = 0
        for line_nm_b_f in nm_b_f.readlines():
            j += 1

            if j < 2:
                continue

            else:
                nm_n_l.append(line_nm_b_f)

    with open(nm_new, "w") as nm_n_f:

        for line_nm_n_l in nm_n_l:
            
            if line_nm_n_l.rstrip():

                nm_n_f.write(line_nm_n_l)

# ==================================================================

# ==================================================================
def nmerbuilder(nmertype, rcut):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter."""

    # Function reused for different types of N-mers.

    if nmertype == "dimers":
        print("Merging monomers with monomers to obtain dimers.\n")
        nmerpatt = "1-*.xyz"
        ucnmlbls = "Dimer"
        numnmlbl = "2"

    elif nmertype == "trimers":
        print("Merging dimers with monomers to obtain trimers.\n")
        nmerpatt = "2-*+*.xyz"
        ucnmlbls = "Trimer"
        numnmlbl = "3"

    elif nmertype == "tetramers":
        print("Merging trimers with monomers to obtain tetramers.\n")
        nmerpatt = "3-*+*+*.xyz"
        ucnmlbls = "Tetramer"
        numnmlbl = "4"
    
    else:
        print_error("The N-mer type must be defined as 'dimers', 'trimers', or 'tetramers'.")

    directory = os.getcwd()
    nmerfs = os.listdir(directory)
    nfrgcntrcll = 12 # TODO: Number of fragments in the central cell for monomer-in-the-cell cutoff.
    moidx = 0
    nmidx = 0
    generatd = 0
    dscrdexs = 0
    dscrdsep = 0

    for monomer in nmerfs:
        
        if (re.match('^1-[0-9]+.xyz$', monomer) and moidx < nfrgcntrcll): # Dimer with at least one monomer in the central cell filter

            for nmer in nmerfs:

                if fnmatch.fnmatch(nmer, nmerpatt):
                    newnm = numnmlbl + "-" + str(monomer)[2:-4] + "+" + str(nmer)[2:]

                    if monomer < nmer: # WARNING: This may not work for trimers, tetramers, ...
                        nmidx += 1
                        
                        mindist = rmin(nmer, monomer)
 
                        if rcut <= mindist: # Separation cutoff filter.
                            print("%s %i (%s) discarded: Separation %3.2f A longer than cutoff %3.2f A" % (ucnmlbls, nmidx, newnm, mindist, rcut))
                            dscrdsep += 1

                            # TODO: Nuclear repulsion energy filter.
                            
                            # TODO: Filter for symetric structures.
 
                        else:
                            nmermerger(monomer, nmer, newnm)
                            print("%s %i (%s) generated: Merged %s and %s" % (ucnmlbls, nmidx, newnm, monomer, nmer))
                            generatd += 1
                    
                    elif nmer == monomer: # Identic structures filter
                        continue

                    else: # Double-counting filter
                        nmidx += 1
                        print("%s %i (%s) discarded: %s of %s and %s already exists" % (ucnmlbls, nmidx, newnm, ucnmlbls, nmer, monomer))
                        dscrdexs += 1

            if (nmertype == "dimers"):
                moidx += 1		

    print("\nTotal number of monomers: %i" % total_monomers)
    print("\nMaximum possible number of dimers: %i" % total_dimers)
    print("\nMaximum possible number of trimers: %i" % total_trimers)
    print("\nMaximum possible number of tetramers: %i" % total_tetramers)
    print("\nDiscarded %i far-separated N-mers and %i to avoid double counting of structures." % (dscrdsep, dscrdexs))
    print("\nGenerated %i N-mers.\n" % generatd)

# ==================================================================

# ==================================================================
def main():
    "Takes a CIF file and computes the crystal lattice energy using a manybody expansion approach."
    
    print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
    print("                             CrystaLattE                             \n")
    print(" The tool for the automated calculation of crystal lattice energies. \n")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

    # ------------------------------------------------------------------
    # Read a CIF file and generate the unit cell.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Run the unit cell through auto_fragments() to extract the number
    # of molecules in the central unit cell.
    # UnitCFrags = auto_fragments(molecule = UnitCell)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Read a CIF file and generates a supercell.
    
    ReadCifIn = 'Benzene-138K.cif'  # CIF input file. # WARNING: Hardcoded for now!
    ReadCifOut = 'Benzene-138K.xyz' # XYZ output file. # WARNING: Hardcoded for now!
    ReadCifA = '5'                  # X replicas. # WARNING: Hardcoded for now!
    ReadCifB = '5'                  # Y replicas. # WARNING: Hardcoded for now!
    ReadCifC = '5'                  # Z replicas. # WARNING: Hardcoded for now!
    
    args = ['', '-i', ReadCifIn, '-o', ReadCifOut, '-b', ReadCifA, ReadCifB, ReadCifC]
    
    print("The following arguments will be passed to the CIF reader script:\n")
    print("./Read_CIF.py" + " ".join(str(arg) for arg in args) + "\n")
    
    Read_CIF.main(args)
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Take the supercell .xyz file
    # And generate .xyz files with all possible monomers.
    print ("")
    
    # Read the lines in the .xyz file recently generated by Read_CIF.py
    with open(ReadCifOut) as fxyz:
        sxyz = fxyz.readlines()
    
    # Generate a SuperCell object.
    #SuperCell = psi4.geometry('\n'.join(sxyz[2:]))
    #SuperCell.update_geometry()
    #print (SuperCell.print_out())
    
    # Generate fragments from SuperCell.
    #CellFrags = auto_fragments(molecule = SuperCell)
    #print (CellFrags.natom())
    #print (CellFrags.nfragments())
    
    # Read the output of the automatic fragmentation.
    
    print("Producing .xyz files for each fragment of the supercell.\n")

    #p4frag = "bztest.p4" # NOTE: for test only!
    p4frag = "bzfrag.p4" # WARNING: Name of the fragmented super cell file is temporarily hardcoded.
    
    numfrags = 664       # WARNING: Number of fragments temporarily hardcoded!
    numfatoms = 12       # WARNING: Number of atoms per fragment hardcoded! What if there are two types of molecules?
    frg_separator = "--" # Fragment separator string.

    psifrg2xyz(p4frag, numfrags, numfatoms, frg_separator)

    # Discard fragments that are not a complete molecule.
    
    print("Detecting fragments with incomplete molecules.\n")

    frags_filter(numfatoms)

    # Organize fragments according to their distance to the center of
    # the supercell, produce molecule files.
    
    print("Organizing molecules according to their separation to the")
    print("center of the supercell.\n")

    frgs2mols()
    
    # Create monomer files from fragment files.

    print("Creating monomers from molecule files.\n")

    mols2mons()
    
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all monomers and generate dimers with all other
    # monomers.
    #
    # Filter dimers that are too distant apart.
    
    rdimer = 3.0
    nmerbuilder("dimers", rdimer)

    # Filter out and keep count of all non-unique dimers, using the
    # nuclear repulsion energy criteria.
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all dimers and generate trimers with all other
    # monomers.
    #
    # Filter trimers that are too distant apart.

    rtrimer = 3.0
    nmerbuilder("trimers", rtrimer)

    # Filter out and keep count of all non-unique trimers, using 
    # ArbAlign.
    # ------------------------------------------------------------------
    
    # .
    # .
    # .
    
    # ------------------------------------------------------------------
    # Run plesantly parallel PSI4 computations on all the final list of 
    # monomers, dimers, trimers, etc.
    #
    # Multiply the resulting energy of each one by the degeneracy factor.
    #
    # Sum results to get a lattice energy.
    # ------------------------------------------------------------------

    print("Execution terminated. Thanks for using CrystaLattE.\n")
# ==================================================================

if __name__ == "__main__":
    main()
