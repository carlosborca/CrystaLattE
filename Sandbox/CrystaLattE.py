#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice 
# energies.
#
# Copyright (c) 2017-2018 
# Carlos H. Borca, Lori A. Burns, and C. David Sherrill.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of CrystaLattE.
#
# CrystaLattE is free software; you can redistribute it and/or modify
# it under the tesms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, version 3.
#
# CrystaLattE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with CrystaLattE; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
# 02110-1301 USA.
#
# @END LICENSE
#

# Import standard Python modules.
import itertools
import multiprocessing
import numpy as np
import os
import sys
import time

# Import outsourced code.
sys.path.insert(0, "Read_CIF")
import Read_CIF

# Import parts of Psi4.
import psi4
from psi4.driver import qcdb
from psi4.driver.qcdb.align import B787
from psi4.driver.qcdb.bfs import BFS
from psi4.driver.qcdb.periodictable import el2mass
from psi4.driver.qcdb.periodictable import el2z


# ======================================================================
def input_parser(in_f_name):
    """.
    """

    keywords = {}
    keywords["nmers_up_to"] = 2 
    keywords["read_cif_output"] = "sc.xyz"
    keywords["read_cif_a"] = 5
    keywords["read_cif_b"] = 5
    keywords["read_cif_c"] = 5
    keywords["r_cut_com"] = 10.0 
    keywords["r_cut_monomer"] = 12.0 
    keywords["r_cut_dimer"] = 10.0 
    keywords["r_cut_trimer"] = 8.0
    keywords["r_cut_tetramer"] = 6.0 
    keywords["r_cut_pentamer"] = 4.0 
    keywords["cle_run_type"] = ["psi4", "quiet"]
    keywords["psi4_method"] = "HF/STO-3G"
    keywords["psi4_memory"] = "500 MB"
    keywords["verbose"] = 1

    with open(in_f_name, "r") as input_file:
    
        for line_inp_f in input_file:
            split_line = line_inp_f.split("=")
            keyword_name = split_line[0].lower()
            keyword_value = split_line[1].split('\n')[0]
       
            if keyword_name == "cle_run_type":
                keyword_value = keyword_value.lower()
                keyword_value = keyword_value.replace(" ", "").split("+")
                
            try:
                keyword_value = float(keyword_value)

            except:
                pass
        
            keywords[keyword_name] = keyword_value

    main(keywords["read_cif_input"], 
         keywords["read_cif_output"],
         keywords["read_cif_a"], 
         keywords["read_cif_b"],
         keywords["read_cif_c"], 
         keywords["nmers_up_to"], 
         keywords["r_cut_com"], 
         keywords["r_cut_monomer"], 
         keywords["r_cut_dimer"], 
         keywords["r_cut_trimer"], 
         keywords["r_cut_tetramer"], 
         keywords["r_cut_pentamer"], 
         keywords["cle_run_type"],
         keywords["psi4_method"],
         keywords["psi4_memory"],
         keywords["verbose"])
# ======================================================================


# ======================================================================
def read_cif_driver(read_cif_input, read_cif_output, read_cif_a, read_cif_b, read_cif_c, verbose=1):
    """Takes strings with the names of a CIF input file and a .xyz 
    output file, as well as integers for the number of replicas of the
    rectangular cell in each direction (A, B, and C). It then calls 
    Read_CIF() and passes that information as arguments to generate an
    .xyz file of the supercell.
    """

    read_cif_arguments = ["", "-i", read_cif_input, "-o", read_cif_output, "-b", read_cif_a, read_cif_b, read_cif_c]

    if verbose >= 2:
        print("\nGenerating the supercell .xyz file.")
        print("\nThe following arguments will be passed to the CIF reader script:")
        print("./Read_CIF.py" + " ".join(str(read_cif_argument) for read_cif_argument in read_cif_arguments) + "\n")
    
    return read_cif_arguments
# ======================================================================


# ======================================================================
def center_supercell(read_cif_output, verbose=0):
    """Takes a string with the filename of the supercell xyz file 
    produced by Read_CIF. Then the center of the supercell coordinates 
    is computed to translate the supercell to the origin.
    
    Returns `scell_geom_max_coords`, which is an np.array with 3 
    numbers: the maximum value of the coordinates x, y, and z of the
    centered supercell.

    It also returns scell_geom and scell_elem, two Numpy arrays with
    geometry and element symbols of the centered supercell, 
    respectively.
    """
    
    if verbose >= 2:
        print("Generating monomers for all complete molecules in the supercell:")

    # Creates two NumPy arrays: one with the coordinates of atoms in the
    # supercell and other with the element symbols of the atoms in it.
    scell_geom = np.loadtxt(read_cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64)
    scell_elem = np.loadtxt(read_cif_output, skiprows=2, usecols=(0), dtype="str")
    
    # Distnaces will be handled in Bohr.
    scell_geom = scell_geom / qcdb.psi_bohr2angstroms

    # Calculation of the supercell center as the midpoint of all
    # coordinates.
    scell_cntr = (np.max(scell_geom, axis=0) - np.min(scell_geom, axis=0))/2.0
    
    if verbose >= 2:
        print("\nCenter of the supercell located at:")
        print("x = %10.5f" % (scell_cntr[0]))
        print("y = %10.5f" % (scell_cntr[1]))
        print("z = %10.5f" % (scell_cntr[2]))
        print("\nThe supercell will be translated to the origin.")

    # Translate the supercell to the origin.
    scell_geom -= scell_cntr
    
    # Return a NumPy array with 3 numbers: the maximum of each x, y, and
    # z.
    scell_geom_max_coords = np.max(scell_geom, axis=0)

    return scell_geom_max_coords, scell_geom, scell_elem
# ======================================================================


# ======================================================================
def supercell2monomers(read_cif_output, r_cut_monomer, verbose=1):
    """Takes a string with the filename of the supercell xyz file
    produced by Read_CIF, and passes it to the `center_supercell()`
    function which translates the supercell to the origin.

    The centered supercell geometries and elements arrays are passed to
    the Breadth-First Search of which returns all fragments found in
    the supercell.

    This function also takes a second argument, a float, with a cutoff
    distance, measured from the origin, which is then used to decide if
    a monomer should be considered or not based on its proximity to the
    origin.

    Returns a dictionary with all the fragments that are within the
    cutoff criterion.
    """
    
    # Centering the supercell.
    scell_geom_max_coords, scell_geom, scell_elem = center_supercell(read_cif_output, verbose)

    # Check if each of dimensions of the supercell satisfies the 
    # condition of being the twice as long as the cutoff. This helps
    # filtering out fragments that contain incomplete molecules located
    # at the edges of the supercell.
    if (r_cut_monomer / qcdb.psi_bohr2angstroms) > np.min(scell_geom_max_coords):
        print("\nWARNING: Cutoff (%3.2f A) longer than half the smallest dimension of the supercell (%3.2f A)." \
              % (r_cut_monomer, np.min(scell_geom_max_coords)*qcdb.psi_bohr2angstroms))
        print("         Please increase the size of the supercell to at least twice r_cut_monomer or reduce the lenght of the cutoff.")

    # Passes the supercell geometry and elements to the breadth-first
    # search algorithm of QCDB to obtain fragments.
    fragments = BFS(scell_geom, scell_elem)
    
    # Two lists containing geometries and elements of a fragment.
    frag_geoms = [scell_geom[fr] for fr in fragments]
    frag_elems = [scell_elem[fr] for fr in fragments]

    # List containing the magnitude of the shortest position vector of
    # each fragment.
    frag_r_mins = []

    for frag in frag_geoms:
        # Get the minium of the norm of the position vectors of all
        # atoms in the fragment.
        r_min_atom = np.min(np.linalg.norm(frag, axis=1))
        frag_r_mins.append(r_min_atom)
    
    # Convert the list of minimums to a NumPy array
    frag_r_mins = np.array(frag_r_mins)

    # Array mapping the indices of the fragments sorted by r_min_atom to
    # the order in the frag_r_mins list.
    mapper = frag_r_mins.argsort()

    # A dictionary for storing N-mers will be created. 

    # Then, for each fragment in the frag_geoms list, an index will be 
    # defined based on the mapper list.
    
    # If the shortest position vector is shorter that the monomers 
    # cutoff, a new dictionary will be created to store all the 
    # information corresponding to such N-mer.
    nmers = {}

    for i in range(len(frag_geoms)):
        index = mapper[i]

        # Discard edges of the supercell and keeps the monomers within a
        # sphere of the cutoff radius.
        if frag_r_mins[index] <= (r_cut_monomer / qcdb.psi_bohr2angstroms):
            name = "1mer-" + str(i)
            
            # The dictionary of one N-mer.
            nmers[name] = {}
            
            nmers[name]["monomers"] = [i]
            nmers[name]["elem"] = frag_elems[index]
            nmers[name]["coords"] = frag_geoms[index]
            nmers[name]["rmin"] = frag_r_mins[index]
            nmers[name]["delimiters"] = []
            nmers[name]["com"] = center_of_mass(frag_elems[index], frag_geoms[index])

    return nmers
# ======================================================================


# ======================================================================
def create_nmer(nmers, ref_monomer, other_monomers, verbose=1):
    """Takes a `nmers` dictionary, and two strings with the keys of a
    refrence monomer and other monomer.

    This function will merge both monomers and create a new N-mer.

    It returns a string `nm_new_name` and a `nm_new` dictionary
    containing the information of the newly generated N-mer.
    """
    
    nm_new = {}

    # Monomers included in the new N-mer.
    nm_new_monomers = []
    nm_new_monomers.append(ref_monomer["monomers"][0])
    nm_new_monomers.extend(other_monomers)
    nm_new["monomers"] = nm_new_monomers
    
    # Elements and coordinates of the atoms in the new N-mer.
    nm_new_elem_arrays = []
    nm_new_coords_arrays = []
    nm_new["atoms_per_monomer"] = []

    for monomer in nm_new_monomers:
        name = "1mer-" + str(monomer)
        nm_new["atoms_per_monomer"].append(len(nmers[name]["elem"]))
        nm_new_elem_arrays.append(nmers[name]["elem"])
        nm_new_coords_arrays.append(nmers[name]["coords"])

    nm_new["elem"] = np.concatenate(nm_new_elem_arrays, axis=0)
    nm_new["coords"] = np.concatenate(nm_new_coords_arrays, axis=0)
    
    # Indices of each monomer in the array of elements and coordinates 
    # (For spliting N-mers into monomers).
    nm_new["delimiters"] = np.cumsum(nm_new["atoms_per_monomer"])
    
    # Norm of the shortest position vector of an atom in the new N-mer.
    # The reference is always closest to the origin so it is always the
    # minimum rmin
    nm_new["rmin"] = ref_monomer["rmin"]

    # Nuclear repulsion energy of the new N-mer.
    nm_new["nre"] = nre(nm_new["elem"], nm_new["coords"])

    # Non additive N-body energy of the N-mer.
    nm_new["nambe"] = 0.0

    # Counter of the number of replicas of this N-mer.
    nm_new["replicas"] = 1

    # Contribution of unique structure to crystal lattice energy.
    nm_new["contrib"] = 0.0

    # Shortest separation vector between the atoms of the monomers in
    # the N-mer.
    nm_new["min_monomer_separations"] = []

    # Separation vector between center of mass of the monomers in the
    # N-mer.
    nm_new["com_monomer_separations"] = []

    for a in nm_new_monomers:
        a_name = "1mer-" + str(a)
        
        for b in nm_new_monomers:
            
            if b <= a:
                continue
            
            b_name = "1mer-" + str(b)

            distm, r_min = distance_matrix(nmers[a_name]["coords"], nmers[b_name]["coords"]) 
            nm_new["min_monomer_separations"].append(r_min)

            nm_new["com_monomer_separations"].append(np.linalg.norm(nmers[a_name]["com"] - nmers[b_name]["com"]))
    
    
    # Criterion to launch energy calculations.
    nm_new["priority"] = 0.0

    priority = 1.0
    
    #for r in nm_new["com_monomer_separations"]:
    for r in nm_new["com_monomer_separations"]:
        one_over_r3 = 1.0/r**3
        priority *= one_over_r3
    
    nm_new["priority"] = priority

    # Key of the new N-mer in the nmers dictionary.
    nm_new_name = str(len(nm_new_monomers)) + "mer-" + "+".join(map(str, nm_new_monomers))

    return nm_new_name, nm_new
# ======================================================================


# ======================================================================
def center_of_mass(elems, geoms):
    """Takes two Numpy arrays, one with the element symbols and one with
    coordinates of a set of atoms and returns a Numpy array with the 
    computed center of mass of the molecule.
    """
    
    com = np.array([0.0, 0.0, 0.0])
    total_mass = 0.0

    for at in range(len(geoms)):
        m = el2mass[elems[at]]
        com += m*geoms[at]
        total_mass += m
    
    com /= total_mass
    
    return com
# ======================================================================


# ======================================================================
def distance_matrix(a, b):
    """Euclidean distance matrix between rows of arrays `a` and `b`.
    Equivalent to `scipy.spatial.distance.cdist(a, b, 'euclidean')`.
    Returns a.shape[0] x b.shape[0] array. It also returns the shortest
    distance between the atoms of `a` and of `b`.
    """

    assert(a.shape[1] == b.shape[1])
    distm = np.zeros([a.shape[0], b.shape[0]])
    
    for i in range(a.shape[0]):
    
        for j in range(b.shape[0]):
            distm[i, j] = np.linalg.norm(a[i] - b[j])

    r_min = np.min(distm)

    return distm, r_min
# ======================================================================


# ======================================================================
def nre(elem, geom):
    """Takes two Numpy arrays, one with the element symbols and one with
    coordinates of a set of atoms and returns a float number with the 
    computed nuclear repulsion energy.
    """
    
    nre = 0.
    for at1 in range(geom.shape[0]):

        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += el2z[elem[at1]] * el2z[elem[at2]] / dist

    return nre
# ======================================================================


# ======================================================================
def build_nmer(nmers, total_monomers, nmer_type, nmer_separation_cutoff, coms_separation_cutoff, verbose=1):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter.
    """

    # Function reused for different types of N-mers.

    if nmer_type == "dimers":
        nm_dictname_pattern = "1mer-"
        num_monomers = 2

    elif nmer_type == "trimers":
        nm_dictname_pattern = "2mer-"
        num_monomers = 3

    elif nmer_type == "tetramers":
        nm_dictname_pattern = "3mer-"
        num_monomers = 4

    elif nmer_type == "pentamers":
        nm_dictname_pattern = "4mer-"
        num_monomers = 5
    
    else:
        print("\nERROR: The N-mer type must be defined as 'dimers', 'trimers', 'tetramers', or 'pentamer'.")

    if verbose >= 2:
        print("")

    counter_new_nmers = 0 # Number of new N-mers generated
    counter_dscrd_sep = 0 # N-mers filtered out by atomic separation
    counter_dscrd_com = 0 # N-mers filtered out by COM separation
    counter_dscrd_rep = 0 # N-mers filtered out as a replicas

    new_nmers = {}

    # TODO: Support for crystals with more than one molecule in the
    #       primitive unit cell.
    num_ref_monomers = 1

    for ref_monomer_idx in range(num_ref_monomers):
        monomer_key = "1mer-" + str(ref_monomer_idx)
        ref_monomer = nmers[monomer_key]

        for other_monomers in itertools.combinations(range(ref_monomer_idx+1, total_monomers), num_monomers-1):

            new_nmer_name, new_nmer = create_nmer(nmers, ref_monomer, other_monomers, verbose)

            max_mon_sep = max(new_nmer["min_monomer_separations"])
            max_com_sep = min(new_nmer["com_monomer_separations"])

            if max_mon_sep > (nmer_separation_cutoff / qcdb.psi_bohr2angstroms):
                
                if verbose >= 2: 
                    print("%s discarded: Atomic separation %3.2f A, longer than cutoff %3.2f A" \
                          % (new_nmer_name, max_mon_sep*qcdb.psi_bohr2angstroms, nmer_separation_cutoff))
                        
                    counter_dscrd_sep += 1
            
            elif max_com_sep > (coms_separation_cutoff / qcdb.psi_bohr2angstroms):
                
                if verbose >= 2:
                    print("%s discarded: Centers of mass separation %3.2f A, longer than cutoff %3.2f A" \
                          % (new_nmer_name, max_com_sep*qcdb.psi_bohr2angstroms, coms_separation_cutoff))

                    counter_dscrd_com += 1

            else:
                found_duplicate = False

                for kexisting, existing in new_nmers.items():

                    # Nuclear repulsion energy filter.
                    if abs(existing["nre"] - new_nmer["nre"]) < 1.e-5:

                        if verbose >= 3:
                            print(kexisting)
                            print(existing["elem"])
                            print(existing["coords"])
                        
                        if verbose >= 3:
                            print(new_nmer_name)
                            print(new_nmer["elem"])
                            print(new_nmer["coords"])
                        
                        # Redirect B787 output.
                        #B787out = kexisting + "-B787.dat"
                        #sys.stdout = open(B787out, 'w')

                        # Block B787 printout
                        sys.stdout = open(os.devnull, 'w')
                        
                        # Block NumPy divide over zero warning printout.
                        np.seterr(divide='ignore')

                        # Call the dreamliner from QCDB.
                        rmsd, mill = B787(rgeom=existing["coords"], cgeom=new_nmer["coords"], runiq=existing["elem"], cuniq=new_nmer["elem"], run_mirror=True, verbose=2)

                        # Reanable printout
                        sys.stdout = sys.__stdout__

                        if rmsd < 1.e-3:
                            found_duplicate = True

                            if verbose >= 2:
                                print("%s discarded: This is a replica of %s." % (new_nmer_name, kexisting))

                            existing["replicas"] += 1
                            counter_dscrd_rep += 1

                            break

                if not found_duplicate:
                    new_nmers[new_nmer_name] = new_nmer

                    if verbose >= 2:
                        print("%s generated: NRE = %.12f au." \
                                % (new_nmer_name, new_nmer["nre"]))

                    counter_new_nmers += 1
    
    if verbose >= 2:
        print("\n{} unique {} were found and generated.".format(counter_new_nmers, nmer_type))

    if verbose >= 2:
        print("\n{} {} did not meet the atomic separation cutoff and were discarded.".format(counter_dscrd_sep, nmer_type))
        print("{} {} did not meet the center of mass separation cutoff and were discarded.".format(counter_dscrd_com, nmer_type))
        print("{} {} were duplicates of another dimer and were discarded.".format(counter_dscrd_rep, nmer_type))
    
    nmers.update(new_nmers)

    return nmers
# ======================================================================


# ======================================================================
def nmer2psiapimol(nmers, knmer, nmer, verbose=0):
    """Takes the `nmers` dictionary; `knmer`, the key of a given N-mer of
    the N-mers dictionary; and `nmer`, its corresponding dictionary.
    Returns a string `psi_api_molecule` that defines a molecule in PSI4
    API mode. This string starts with three lines specifying the use of
    atomic units, avoidance of center of mass translation, and
    avoidance of reorientation. The next lines contain the element
    symbol and coordinates, for each atom in the N-mer, separating each
    monomer with lines containing a double hyphen.
    """

    psi_api_molecule = "\nunits = au\n"
    psi_api_molecule += "no_com\n"
    psi_api_molecule += "no_reorient\n"

    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psi_api_molecule += "--\n"

        psi_api_molecule += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    return psi_api_molecule
# ======================================================================


# ======================================================================
def nmer2psithon(nmers, knmer, nmer, psi4_method, psi4_memory, verbose=0):
    """.
    """
    
    psithon_input =  "# Psi4 Psithon input file for N-mer: {}\n".format(knmer)
    
    psithon_input += "\nmemory {}\n".format(psi4_memory)
    
    mymol = knmer.replace("2mer", "Dimer").replace("3mer", "Trimer").replace("4mer", "Tetramer").replace("5mer", "Pentamer").replace("-", "_").replace("+", "_")
    psithon_input += "\nmolecule {} {{\n".format(mymol)
    
    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psithon_input += "--\n"
        
        psithon_input += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    psithon_input += "units = au\n"
    psithon_input += "no_com\n"
    psithon_input += "no_reorient\n"
    
    psithon_input += "}\n"

    psithon_input += "\nset {\n"
    psithon_input += "  e_convergence = 8\n"
    psithon_input += "  scf_type df\n"
    psithon_input += "  mp2_type df\n" 
    psithon_input += "  cc_type df\n"
    psithon_input += "  freeze_core true\n"
    psithon_input += "}\n"

    # Hartree-Fock is called with the 'scf' string in Psithon mode.
    if psi4_method.lower().startswith("hf"):
        psithon_method = "scf" + psi4_method[2:]

    else:
        psithon_method = psi4_method

    # This is a temporary fix because the N-Body wrapper executed 4 calculations for a dimer when doing VMFC!!
    if len(nmer["monomers"]) > 2:
        psithon_input += "\nenergy('{}', bsse_type = 'vmfc')\n".format(psithon_method)

    else:
        psithon_input += "\nenergy('{}', bsse_type = 'cp')\n".format(psithon_method)

    print("\nPSI4 Molecule of %s:\n" % knmer)
    print(psithon_input)

    return psithon_input
# ======================================================================


# ======================================================================
def psi4_energies(nmers, knmer, nmer, cpus, cle_run_type, psi4_method, psi4_memory, verbose=0):
    """.
    """

    # If the output is going to be kept, setup the filename.
    if "quiet" in cle_run_type:
        psi4.core.be_quiet()
        
    # If the output is not kept, do not print to screen.
    else:
        p4out = knmer + ".dat"
        psi4.core.set_output_file(p4out)

    psi_api_molecule = nmer2psiapimol(nmers, knmer, nmer, verbose)
    mymol = psi4.geometry(psi_api_molecule)
    
    # Set the number of threads to run Psi4.
    psi4.core.set_num_threads(cpus)
    
    # Set Psi4 memory.
    psi4.set_memory(psi4_memory)
    
    # Set Psi4 block of options.
    psi4.set_options({'scf_type': 'df', 'mp2_type': 'df', 'cc_type': 'df', 'freeze_core': 'true', 'e_convergence': '8'})
    
    # Execute Psi4 energy calculations, unless running on test mode.
    
    # Example:  psi4.energy('MP2/aug-cc-pV[D,T]Z', molecule=he_tetramer, bsse_type=['cp', 'nocp', 'vmfc'])
    #           psi4.energy('HF/STO-3G', molecule=mymol, bsse_type=['vmfc'])
    #           psi4.energy('MP2/aug-cc-pVDZ', molecule=mymol, bsse_type=['vmfc'])
    
    if "test" not in cle_run_type:
        
        # This is a temporary fix because the N-Body wrapper executed 4 calculations for a dimer when doing VMFC!!
        if len(nmer["monomers"]) > 2:
            psi4.energy(psi4_method, molecule=mymol, bsse_type=['vmfc'])
        
        else:
            # psi4.energy(psi4_method, molecule=mymol, bsse_type=['vmfc'])
            psi4.energy(psi4_method, molecule=mymol, bsse_type=['cp'])
    
    # Get the non-additive n-body contribution, exclusive of all
    # previous-body interactions.
    varstring = "CP-CORRECTED " + str(len(nmer["monomers"])) + "-BODY INTERACTION ENERGY"
    
    # This comment and the line above is a temporary fix because the N-Body wrapper executed 4 calculations for a dimer!!
    # varstring = "VMFC-CORRECTED " + str(len(nmer["monomers"])) + "-BODY INTERACTION ENERGY"
    
    n_body_energy = psi4.core.get_variable(varstring)
    
    if len(nmer["monomers"]) > 2:
        varstring = "VMFC-CORRECTED " + str(len(nmer["monomers"])-1) + "-BODY INTERACTION ENERGY"
        n_minus_1_body_energy = psi4.core.get_variable(varstring)
        nmer["nambe"] = n_body_energy - n_minus_1_body_energy

    else:
        nmer["nambe"] = n_body_energy
# ======================================================================


# ======================================================================
def cle_manager(nmers, cle_run_type, psi4_method, psi4_memory, verbose=0):
    """.
    """

    global crystal_lattice_energy
    global results

    crystal_lattice_energy = 0.0
    results = []
    
    for knmer, nmer in nmers.items():

        # Energies are not calculated for monomers. Rigid body approximation.
        if len(nmer["monomers"]) == 1:
            continue

        # Find out the number of CPUs in the local system.
        cpus = multiprocessing.cpu_count()
        
        # Start wall-clock timer.
        energies_start = time.time()
            
        # Run energies in PSI4.
        if "psi4" in cle_run_type:
            psi4_energies(nmers, knmer, nmer, cpus, cle_run_type, psi4_method, psi4_memory, verbose)

        # Stop wall-clock timer.
        energies_end = time.time()
        
        # Calculate execution time.
        energies_wallclock = energies_end - energies_start
    
        nmer["contrib"] = nmer["nambe"] * nmer["replicas"] / float(len(nmer["monomers"]))
        crystal_lattice_energy += nmer["contrib"]
        
        rcomseps = ""

        for r in nmer["min_monomer_separations"]:
            rcomseps += "{:6.3f} ".format(r * qcdb.psi_bohr2angstroms)

        nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                knmer, 
                nmer["nambe"] * qcdb.psi_hartree2kcalmol * qcdb.psi_cal2J, 
                nmer["replicas"], 
                nmer["contrib"] * qcdb.psi_hartree2kcalmol * qcdb.psi_cal2J,
                crystal_lattice_energy * qcdb.psi_hartree2kcalmol * qcdb.psi_cal2J,
                nmer["priority"],
                rcomseps)
        
        results.append(nmer_result)

        if verbose >= 2:
            print("{} elapsed {:.2f} seconds processing on {} threads. Additive Lattice Energy = {:9.8f} KJ/mol".format(
                knmer, 
                energies_wallclock, 
                cpus, 
                crystal_lattice_energy * qcdb.psi_hartree2kcalmol * qcdb.psi_cal2J))
        
    return crystal_lattice_energy
# ======================================================================


# ======================================================================
def print_results(verbose=0):
    """.
    """

    if verbose >= 1:
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+-----------------")
        print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (KJ/mol) |  (#) |     (KJ/mol) |      (KJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+-----------------")
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+-----------------")
        print("\nCrystal Lattice Energy (Eh)       = {:5.8f}".format(crystal_lattice_energy))
        print("Crystal Lattice Energy (KJ/mol)   = {:9.8f}".format(crystal_lattice_energy * qcdb.psi_hartree2kcalmol * qcdb.psi_cal2J))
        print("Crystal Lattice Energy (Kcal/mol) = {:9.8f}\n".format(crystal_lattice_energy * qcdb.psi_hartree2kcalmol))
# ======================================================================


# ======================================================================
def print_end_msg(start, verbose=0):
    """.
    """
    
    if verbose >= 1:
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("Execution terminated succesfully.")
        print("Total elapsed wall-clock time: {:.2f} seconds\n".format(time.time() - start))
        print("Thank you for using CrystaLatte.\n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

# ======================================================================


# ======================================================================
def main(read_cif_input, read_cif_output="sc.xyz", read_cif_a=5, read_cif_b=5, read_cif_c=5, nmers_up_to=2, r_cut_com=10.0, r_cut_monomer=12.0, r_cut_dimer=10.0, r_cut_trimer=8.0, r_cut_tetramer=6.0, r_cut_pentamer=4.0, cle_run_type=["psi4", "quiet"], psi4_method="HF/STO-3G", psi4_memory="500 MB", verbose=1):
    """Takes a CIF file and computes the crystal lattice energy using a
    many-body expansion approach.
    """
   
    # Start counting time.
    start = time.time()
    
    # Check proper input filename.
    if read_cif_input.endswith(".cif"):
        outf = read_cif_input[:-4] + ".out"
    
    else:
        print("CrystaLattE needs a .cif file as input file.")

    # Print program header.
    if verbose >= 1:
        print("")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("                              CrystaLattE                              \n")
        print("  The tool for the automated calculation of crystal lattice energies.  \n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

    # Read a CIF file and generate the unit cell.
    read_cif_arguments = read_cif_driver(read_cif_input, read_cif_output, read_cif_a, read_cif_b, read_cif_c, verbose)
    Read_CIF.main(read_cif_arguments)
    
    # Read the output of the automatic fragmentation.
    nmers = supercell2monomers(read_cif_output, r_cut_monomer, verbose)
    total_monomers = len(nmers)

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
        build_nmer(nmers, total_monomers, "dimers", r_cut_dimer, r_cut_com, verbose)

    if nmers_up_to >= 3:
        
        if verbose >= 2:
            print("\nMerging dimers with monomers to obtain trimers.")
        build_nmer(nmers, total_monomers, "trimers", r_cut_trimer, r_cut_com, verbose)

    if nmers_up_to >= 4:
        
        if verbose >= 2:
            print("\nMerging trimers with monomers to obtain tetramers.")
        build_nmer(nmers, total_monomers, "tetramers", r_cut_tetramer, r_cut_com, verbose)

    if nmers_up_to == 5:
        
        if verbose >= 2:
            print("\nMerging tetramers with monomers to obtain pentamers.")
        build_nmer(nmers, total_monomers, "pentamers", r_cut_pentamer, r_cut_com, verbose)

    if nmers_up_to > 5:
        print("\nERROR: The current implementation of CrystaLattE is limited to pentamers.")
        print("       Please use 2 <= nmer_up_to < 5.")
   
    if verbose >= 2:
        print ("\nComputing interaction energies of N-mers:")

    cle_manager(nmers, cle_run_type, psi4_method, psi4_memory, verbose)
    # ------------------------------------------------------------------
    
    # Print the final results.
    print("")
    print_results(verbose)
    
    # Print exit message and timings information.
    print_end_msg(start, verbose)
# ======================================================================


if __name__ == "__main__":

    # Hard-coded Test
    if "CrystaLattE" in sys.argv[-1]:
        main(   read_cif_input="../Tests/Benzene/Benzene.cif",
                read_cif_output="../Tests/Benzene/Benzene.xyz",
                read_cif_a=2,
                read_cif_b=2,
                read_cif_c=2,
                nmers_up_to=2,
                r_cut_com=5.1,
                r_cut_monomer=3.3,
                r_cut_dimer=2.7,
                r_cut_trimer=2.7,
                r_cut_tetramer=2.7,
                r_cut_pentamer=5.6,
                cle_run_type=["psi4", "quiet"],
                psi4_method="HF/STO-3G",
                psi4_memory="500 MB",
                verbose=2)

        print("><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        print("WARNING: No imput was provided. The previous execution was just a test.")
        print("><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")

    # Normal execution using an input file.
    else:
        input_parser(sys.argv[-1])
