#!/usr/bin/env python

import fnmatch
import math
import numpy as np
import re
import os
import sys

# Imports of outsourced code.
sys.path.insert(0, "Read_CIF")
import Read_CIF

# Import parts of Psi4.
import psi4
from psi4.driver import qcdb
from psi4.driver.qcdb.bfs import BFS
from psi4.driver.qcdb.periodictable import el2z
from psi4.driver.qcdb.align import B787

# ==================================================================
def read_cif_driver(read_cif_input, read_cif_output, read_cif_a, read_cif_b, read_cif_c, verbose=1):
    """Takes a the names of a CIF input file, a .xyz output file, and the number of replicas of the rectangular cell in each direction (A, B, and C).
    It then calls Read_CIF() and passes that information as arguments to generate an .xyz file of the supercell."""

    read_cif_arguments = ["", "-i", read_cif_input, "-o", read_cif_output, "-b", read_cif_a, read_cif_b, read_cif_c]

    if verbose >= 2:
        print("\nGenerating the supercell .xyz file.")
        print("\nThe following arguments will be passed to the CIF reader script:")
        print("./Read_CIF.py" + " ".join(str(read_cif_argument) for read_cif_argument in read_cif_arguments))
    
    return read_cif_arguments
# ==================================================================

# ==================================================================
def supercell2monomers(read_cif_output, r_cut_monomer, verbose=1):
    """Take the supercell xyz file produced by Read_CIF and break it into fragments."""
    
    if verbose >= 1:
        print("Generating monomers for all complete molecules in the supercell:")

    scell_geom = np.loadtxt(read_cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64) # Creates a NumPy array with the coordinates of atoms in the supercell
    scell_elem = np.loadtxt(read_cif_output, skiprows=2, usecols=(0), dtype="str") # Creates a NumPy array with the element symbols of the atoms in the supercell
    
    scell_geom = scell_geom / qcdb.psi_bohr2angstroms # Coordinates will be handled in Bohr

    scell_cntr = (np.max(scell_geom, axis=0) - np.min(scell_geom, axis=0))/2.0 # Determines the supercell center
    
    if verbose >= 2:
        print("\nCenter of the supercell located at:")
        print("x = %10.5f" % (scell_cntr[0]))
        print("y = %10.5f" % (scell_cntr[1]))
        print("z = %10.5f" % (scell_cntr[2]))
        print("\nThe supercell will be translated to the origin")

    scell_geom -= scell_cntr # Translates the supercell to the origin

    scell_geom_max_coords = np.max(scell_geom, axis=0) # Returns an np.array with 3 numbers: the maximum of x, of y, and of z

    if verbose >= 2:
        print("\nChecking if the size of the supercell is acceptable")

    if (r_cut_monomer / qcdb.psi_bohr2angstroms) > np.min(scell_geom_max_coords):
        print("\nERROR: Cutoff (%3.2f A) longer than half the smallest dimension of the supercell (%3.2f A)." \
              % (r_cut_monomer, np.min(scell_geom_max_coords)*qcdb.psi_bohr2angstroms))
        print("       Please increase the size of the supercell to at least twice r_cut_monomer or reduce the lenght of the cutoff.")

    fragments = BFS(scell_geom, scell_elem) # Passes the supercell geometry and elements the breadth-first search algorithm of QCDB to obtain fragments
    frag_geoms = [scell_geom[fr] for fr in fragments]
    frag_elems = [scell_elem[fr] for fr in fragments]

    frag_r_mins = [] # List containing the magnitude of the shortest position vector of each fragment.

    for frag in frag_geoms:
        r_min_atom = np.min(np.linalg.norm(frag, axis=1)) # Get the minium of the norm of the position vectors of all atoms in the fragment.
        frag_r_mins.append(r_min_atom)
    
    frag_r_mins = np.array(frag_r_mins) # Conversion to NumPy array.

    mapper = frag_r_mins.argsort() # Array mapping the indices of the fragments sorted by r_min_atom to the order in the frag_r_mins list.

    # A dictionary for storing N-mers will be created.
    # Then, for each fragment in the frag_geoms list, an index will be defined based on the mapper list.
    # If the shortest position vector is shorter that the monomers cutoff,
    # a new dictionary will be created to store all the information corresponding to such N-mer.

    nmers = {} # The N-Mers dictionary

    for i in range(len(frag_geoms)):
        index = mapper[i]

        if frag_r_mins[index] <= (r_cut_monomer / qcdb.psi_bohr2angstroms): # Discard edges of the supercell and keeps the monomers within a sphere of the cutoff radius.
            name = "1mer-" + str(i)
            
            nmers[name] = {} # The dictionary of one N-mer.
            
            nmers[name]["monomers"] = [i]
            nmers[name]["elem"] = frag_elems[index]
            nmers[name]["coords"] = frag_geoms[index]
            nmers[name]["rmin"] = frag_r_mins[index]

            if verbose >= 3:
                print(name)
                print(nmers[name])
    
    return nmers
# ==================================================================
def merge_nmers(nmers, nm_a, nm_b, verbose=1):
    
    nm_new = {}

    # Monomers included in the new N-mer.
    nm_new_monomers = nm_a["monomers"] + nm_b["monomers"]
    nm_new_monomers.sort()
    nm_new["monomers"] = nm_new_monomers
    
    # Elements and coordinates of the atoms in the new N-mer.
    nm_new_elem_arrays = []
    nm_new_coords_arrays = []
    for monomer in nm_new_monomers:
        name = "1mer-" + str(monomer)
        nm_new_elem_arrays.append(nmers[name]["elem"])
        nm_new_coords_arrays.append(nmers[name]["coords"])

    nm_new["elem"] = np.concatenate(nm_new_elem_arrays, axis=0)
    nm_new["coords"] = np.concatenate(nm_new_coords_arrays, axis=0)

    # Norm of the shortest position vector of an atom in the new N-mer.
    nm_new_rmin = min([nm_a["rmin"], nm_b["rmin"]])
    nm_new["rmin"] = nm_new_rmin

    nm_new_name = str(len(nm_new_monomers)) + "mer-" + "+".join(map(str, nm_new_monomers))

    nm_new["nre"] = nre(nm_new["elem"], nm_new["coords"])
    nm_new["replicas"] = 1

    return nm_new_name, nm_new
# ==================================================================

# ==================================================================
def distance_matrix(a, b):
    """Euclidean distance matrix between rows of arrays `a` and `b`. Equivalent to
    `scipy.spatial.distance.cdist(a, b, 'euclidean')`. Returns a.shape[0] x b.shape[0] array."""

    assert(a.shape[1] == b.shape[1])
    distm = np.zeros([a.shape[0], b.shape[0]])
    
    for i in range(a.shape[0]):
    
        for j in range(b.shape[0]):
            distm[i, j] = np.linalg.norm(a[i] - b[j])

    r_min = np.min(distm)

    return distm, r_min
# ==================================================================

# ==================================================================
def nre(elem, geom):
    """Nuclear repulsion energy"""
    
    nre = 0.
    for at1 in range(geom.shape[0]):

        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += el2z[elem[at1]] * el2z[elem[at2]] / dist

    return nre
# ==================================================================

# ==================================================================
def build_nmer(nmers, nmer_type, nmer_separation_cutoff, verbose=1):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter."""

    # Function reused for different types of N-mers.

    if nmer_type == "dimers":
        nm_dictname_pattern = "1mer-"
        nm_txt_lbl = "Dimer"
        nm_num_lbl = "2"

    elif nmer_type == "trimers":
        nm_dictname_pattern = "2mer-"
        nm_txt_lbl = "Trimer"
        nm_num_lbl = "3"

    elif nmer_type == "tetramers":
        nm_dictname_pattern = "3mer-"
        nm_txt_lbl = "Tetramer"
        nm_num_lbl = "4"

    elif nmer_type == "pentamers":
        nm_dictname_pattern = "4mer-"
        nm_txt_lbl = "Pentamer"
        nm_num_lbl = "5"
    
    else:
        print("\nERROR: The N-mer type must be defined as 'dimers', 'trimers', 'tetramers', or 'pentamer'.")

    counter_new_nmers = 0
    counter_dscrd_spi = 0
    counter_dscrd_exs = 0
    counter_dscrd_sep = 0

    new_nmers = {}

    for monomer_key, monomer in nmers.items():

        if monomer_key.startswith("1mer-"):

            for nmer_key, nmer in nmers.items():

                if nmer_key.startswith(nm_dictname_pattern):
                    if (monomer["monomers"][0] in nmer["monomers"]):
                        continue
                    new_nmer_name, new_nmer = merge_nmers(nmers, monomer, nmer, verbose)
                    distm, r_min = distance_matrix(monomer["coords"], nmer["coords"])

                    if r_min <= 1.e-5: # Superimposed monomers filter.
                        if verbose >= 2:
                            print("%s %s discarded: Separation %3.2f A, superimposed structures" \
                                  % (nm_txt_lbl, new_nmer_name, r_min))
                        
                        counter_dscrd_spi += 1

                    elif r_min > (nmer_separation_cutoff / qcdb.psi_bohr2angstroms): # Separation cutoff filter.
                        if verbose >= 2:
                            print("%s %s discarded: Separation %3.2f A, longer than cutoff %3.2f A" \
                                  % (nm_txt_lbl, new_nmer_name, r_min*qcdb.psi_bohr2angstroms, nmer_separation_cutoff))
                        
                        counter_dscrd_sep += 1

                    # TODO: Nuclear repulsion energy filter (?)
                    else:
                        found_duplicate = False
                        for existing in new_nmers.values():
                            
                            if abs(existing["nre"] - new_nmer["nre"]) < 1.e-5:
                                # check if this is superimposable on the other structure or not
                                # if superimposable then discard.
                                
                                print(existing["coords"], new_nmer["coords"], existing["elem"], new_nmer["elem"])
                                
                                # Call the dream(a)li(g)ner of QCDB.
                                rmsd, mill = B787(rgeom=existing["coords"], cgeom=new_nmer["coords"], runiq=existing["elem"], cuniq=new_nmer["elem"])
                                if rmsd < 1.e-3:
                                    found_duplicate = True
                                    print("%s %s discarded: A duplicate of this N-mer was generated before" % (nm_txt_lbl, new_nmer_name))
                                    existing["replicas"] += 1
                                    break
                            else:
                                print("This one does not have matching NRE ", existing["nre"], new_nmer["nre"])

                        if not found_duplicate:
                            new_nmers[new_nmer_name] = new_nmer
                            
                            if verbose >= 2:
                                print("%s %s generated" \
                                      % (nm_txt_lbl, new_nmer_name))
                        
                            counter_new_nmers += 1
    
    nmers.update(new_nmers)

    return nmers
# ==================================================================

# ==================================================================
#def main(input_cif_file, r_cut_monomer=10.0, r_cut_dimer=10.0, r_cut_trimer=10.0, r_cut_tetramer=10.0, r_cut_pentamer=10.0, verbose=1):
def main(read_cif_input, read_cif_output, read_cif_a, read_cif_b, read_cif_c, fragmented_supecell_file, num_frags, num_atoms_frag, fragment_separator, nmers_up_to=3, r_cut_monomer=10.0, r_cut_dimer=10.0, r_cut_trimer=10.0, r_cut_tetramer=10.0, r_cut_pentamer=10.0, verbose=1):
    "Takes a CIF file and computes the crystal lattice energy using a manybody expansion approach."
    
    # Print program header.
    if verbose >= 1:
        print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("                             CrystaLattE                             \n")
        print(" The tool for the automated calculation of crystal lattice energies. \n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

    # Read a CIF file and generate the unit cell.
    read_cif_arguments = read_cif_driver(read_cif_input, read_cif_output, read_cif_a, read_cif_b, read_cif_c, verbose)
    Read_CIF.main(read_cif_arguments)
    
    # Read the output of the automatic fragmentation.
    nmers = supercell2monomers(read_cif_output, r_cut_monomer, verbose)
    
    # Loop through all monomers and monomers in the central unit cell
    # to generate dimers with at least one monomer in the central cell.
    # Then loop through all existing N-mers and generate higher-order 
    # N-mers with all monomers.

    if nmers_up_to < 2:
        print("\nERROR: CrystaLattE is designed to use at least dimers.")
        print("       Please use 2 <= nmer_up_to < 5.")
    
    if nmers_up_to >= 2:
        
        if verbose >= 2:
            print("\nMerging monomers with monomers to obtain dimers.")
        
        build_nmer(nmers, "dimers", r_cut_dimer, verbose)

    if nmers_up_to >= 3:
        
        if verbose >= 2:
            print("\nMerging dimers with monomers to obtain trimers.")

        build_nmer(nmers, "trimers", r_cut_trimer, verbose)

    if nmers_up_to >= 4:
        
        if verbose >= 2:
            print("\nMerging trimers with monomers to obtain tetramers.")

        build_nmer(nmers, "tetramers", r_cut_tetramer, verbose)

    if nmers_up_to == 5:
        
        if verbose >= 2:
            print("\nMerging tetramers with monomers to obtain pentamers.")

        build_nmer(nmers, "pentamers", r_cut_pentamer, verbose)

    if nmers_up_to > 5:
        print("\nERROR: The current implementation of CrystaLattE is limited to pentamers.")
        print("       Please use 2 <= nmer_up_to < 5.")
    
    # ------------------------------------------------------------------
    # TODO: Run plesantly parallel PSI4 computations on all the final list of 
    # monomers, dimers, trimers, etc.
    #
    # TODO: Multiply the resulting energy of each one by the degeneracy factor.
    #
    # TODO: Sum results to get a lattice energy.
    # ------------------------------------------------------------------

# ==================================================================

if __name__ == "__main__":
    #main(read_cif_input, nmers_up_to=3)
    
    # Test up to trimers, with reduced supercell of benzene, and shortened cutoffs.
    main(   read_cif_input="Benzene-138K.cif", 
            read_cif_output="Benzene-138K.xyz", 
            read_cif_a=3, 
            read_cif_b=3, 
            read_cif_c=3, 
            fragmented_supecell_file="bztest.p4", 
            num_frags=664, 
            num_atoms_frag=12, 
            fragment_separator="--", 
            nmers_up_to=2, 
            r_cut_monomer=5.0, 
            r_cut_dimer=3.0, 
            r_cut_trimer=3.0, 
            r_cut_tetramer=3.0, 
            r_cut_pentamer=3.0, 
            verbose=2)

    # Test up to trimers, with water supercell, and shortened cutoffs.
    #main(   read_cif_input="ice-Ih.cif", \
    #        read_cif_output="ice-Ih.xyz", \
    #        read_cif_a=5, \
    #        read_cif_b=5, \
    #        read_cif_c=5, \
    #        fragmented_supecell_file="icefrag.p4", \
    #        num_frags=13, \
    #        num_atoms_frag=3, \
    #        fragment_separator="--", \
    #        nmers_up_to=3, \
    #        r_cut_dimer=3.0, \
    #        r_cut_trimer=3.0, \
    #        r_cut_tetramer=3.0, \
    #        r_cut_pentamer=3.0, \
    #        verbose=2)
