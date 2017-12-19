#!/usr/bin/env python

import fnmatch
import math
import re
import os
import sys

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
def psifragments2xyzs(fragmented_supecell_file, num_frags, num_atoms_frag, fragment_separator):
    """Takes the output of the automatic fragmentation procedure, the total number of fragments,
    the number of atoms in one fragment, and the string that separates the fragments.
    Returns .xyz files for each fragment."""

    # Read the output of the automatic fragmentation.
    with open(fragmented_supecell_file) as supercell_fragments:
        sc_frags_f_lines = supercell_fragments.readlines()

    # Generate .xyz files for each fragment.
    counter_fragments = 1         # Counter of processed fragments.
    counter_lines     = 0         # Counter of lines.
    header_line       = True      # Flag for first line of a fragment .xyz file.

    for sc_frags_f_line in sc_frags_f_lines:
        counter_lines += 1

        if sc_frags_f_line.startswith(fragment_separator):
            counter_fragments += 1
            header_line = True

        else:
            fragment_filename = "f" + str(counter_fragments).zfill(len(str(num_frags))) \
                     + ".xyz"
            with open(fragment_filename, "a") as fragment_xyz_file: # WARNING: File exists?

                if header_line == True:
                    fragment_xyz_file.write(str(num_atoms_frag) + "\n" + "Fragment " + str(counter_fragments) + "\n")
                    header_line = False

                if header_line == False:
                    fragment_xyz_file.write(sc_frags_f_line)
    
    print("%i fragment .xyz files were created." % counter_fragments)

    return
# ==================================================================

# ==================================================================
def fragment_filter(number_atoms_per_frag): # WARNING: The number of atoms in a fragment probably should be passed as a list.
    """It takes the possible numbers of atoms in a fragment and discards all fragments files that contain a smaller number of atoms."""

    # TODO: What if there are more than two types of molecules in the unit cell? 
    # The number of atoms per framgent should be passed as a list with possible numbers of atoms in a fragment.

    directory = os.getcwd()
    files = os.listdir(directory)
    counter_removed_fragments = 0

    for file in files:

        if re.match('^f[0-9]+.xyz$', file): # Match filenames f???.xyz

            with open(file, 'r') as fragment_file:
                frag_file_lines = fragment_file.readlines()

                if len(frag_file_lines) != number_atoms_per_frag + 2: # Match files with unexpected atoms.
                    print("Expected %i atoms. Found only %i. Removing: %s" % (number_atoms_per_frag, (len(frag_file_lines) - 2), file))
                    os.remove(os.path.join(directory,file))
                    counter_removed_fragments += 1

    print("\nRemoved %s fragments.\n" % counter_removed_fragments)

    return
# ==================================================================

# ==================================================================
def supercell_center(frag_filename_pattern):
    """Compute the center of coordinates of .xyz files with a matching pattern in the filename.
    Returns a list with fragment names, x, y, and z coordinates in lists of floats and the coordinates of the center."""

    frag_filenames_list = [] # Fragments file names.
    frag_coords_x = []           # Values of x coordinates.
    frag_coords_y = []           # Values of y coordinates.
    frag_coords_z = []           # Values of z coordinates.

    directory = os.getcwd()
    files = os.listdir(directory)

    for file in files:

        if re.match(frag_filename_pattern, file):
            
            with open(file, 'r') as frag_file:
                frag_file_lines = frag_file.readlines()

                i = 0
                for i in range(len(frag_file_lines)):

                    if i < 2:
                        continue

                    else:
                        frag_info = frag_file_lines[i].split()
                        frag_filenames_list.append(file)
                        frag_coords_x.append(float(frag_info[1]))
                        frag_coords_y.append(float(frag_info[2]))
                        frag_coords_z.append(float(frag_info[3]))

    cntr_x = (max(frag_coords_x) - min(frag_coords_x))/2.0
    cntr_y = (max(frag_coords_y) - min(frag_coords_y))/2.0
    cntr_z = (max(frag_coords_z) - min(frag_coords_z))/2.0

    print("Center of the supercell located at:")
    print("x = %10.5f"   % (cntr_x))
    print("y = %10.5f"   % (cntr_y))
    print("z = %10.5f\n" % (cntr_z))

    return frag_filenames_list, frag_coords_x, frag_coords_y, frag_coords_z, cntr_x, cntr_y, cntr_z 
# ==================================================================

# ==================================================================
def position_vec_mag(frag_filename_pattern):
    """Using the supercell_center() function it computes the magnitude of the position vector of the atoms read by the supercell_center() function.
    Returns the filename and the magnitude of the position vector."""

    frag_fnames_list, frag_x, frag_y, frag_z, cntr_x, cntr_y, cntr_z = supercell_center(frag_filename_pattern)

    atoms_r = []

    for j in range(len(frag_x)):
        atoms_r.append(math.sqrt((frag_x[j] - cntr_x)**2.0 + (frag_y[j] - cntr_y)**2.0 + (frag_z[j] - cntr_z)**2.0))

    print("Magnitudes of the position vectors calculated.\n")

    return frag_fnames_list, atoms_r
# ==================================================================

# ==================================================================
def frag_atoms_r_list(frag_filename_pattern):
    """Using position vector magnitude, computed in position_vec_mag(). This function returns a list of fragments filenames
    with fragments containing atoms closest to the center of the supercell first."""

    frag_filenames_list, r = position_vec_mag(frag_filename_pattern)

    frag_proximity_list = []

    for pair in sorted(enumerate(r), key=lambda item:item[1]):
        if frag_filenames_list[pair[0]] in frag_proximity_list:
            continue
        else:
            frag_proximity_list.append(frag_filenames_list[pair[0]])

    print("Molecules organized based on the closest atom to the center of the supercell.\n")

    return frag_proximity_list
# ==================================================================

# ==================================================================
def fragments2molecules():
    """Takes fragents and a list of proximity to the center of the supercell, generated by the frag_atoms_r_list() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell."""

    directory = os.getcwd()
    frag_prox_ordered_list = frag_atoms_r_list("^f[0-9]+.xyz$")
    counter_new_molecule = 0

    for file in frag_prox_ordered_list:
        counter_new_molecule += 1
        filename_new_molecule = "m" + str(counter_new_molecule).zfill(len(frag_prox_ordered_list[-1]) - 5) + ".xyz"

        with open(file, 'r') as frag_file, open(filename_new_molecule, 'w') as mol_file:

            for frag_file_line in frag_file.readlines():

                if frag_file_line.startswith("Fragment"):
                    frag_file_new_line = "Molecule " + str(counter_new_molecule) + " (" + frag_file_line[:-1] + ")\n"
                    mol_file.write(frag_file_new_line)

                else:
                    mol_file.write(frag_file_line)

        os.remove(os.path.join(directory, file))

        print("Fragment %s has been rewritten to: %s" % (file, filename_new_molecule))

    print("")
# ==================================================================

# ==================================================================
def molecules2monomers():
    """Takes molecules and filters and a list of proximity to the center of the supercell, generated by the frag_atoms_r_list() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell.
    It returns four global variables with the number of monomeers and the maximum possible number of dimers, trimers, tetramers."""

    directory = os.getcwd()
    files = os.listdir(directory)
    monidx = 0
    
    global counter_total_monomers
    counter_total_monomers = 0

    for file in files:
        
        if re.match('^m[0-9]+.xyz$', file):
            monidx += 1
            monfname = "1-" + str(file)[1:]

            with open(file, 'r') as moleculef, open(monfname, 'w') as monomerf:
                counter_total_monomers += 1

                for l in moleculef.readlines():

                    if l.startswith("Molecule"):
                        newl = "Monomer " + str(monidx) + "\n"
                        monomerf.write(newl)

                    else:
                        monomerf.write(l)

            os.remove(os.path.join(directory, file))

            print("Molecule " + str(file) + " has been rewritten to: " + str(monfname))

    print("\nTotal number of mononmers: %i" % counter_total_monomers)
    
    global counter_total_dimers
    counter_total_dimers = (counter_total_monomers)*(counter_total_monomers - 1)/2

    global counter_total_trimers
    counter_total_trimers = (counter_total_dimers)*(counter_total_monomers - 2)/3

    global counter_total_tetramers
    counter_total_tetramers = (counter_total_trimers)*(counter_total_monomers - 3)/4

    global counter_total_pentamers
    counter_total_pentamers = (counter_total_tetramers)*(counter_total_monomers - 4)/5

# ==================================================================

# ==================================================================
def r_min(nmer1, nmer2):
    """Takes two N-mer files and returns separation between closest atoms belonging to different N-mers."""

    with open(nmer1, 'r') as nmer1_file, open(nmer2, 'r') as nmer2_file:
        nmer1_file_lines = nmer1_file.readlines()
        nmer2_file_lines = nmer2_file.readlines()

        nmer1_info_list = []
        nmer2_info_list = []

        atomic_separations = []

        i = 0
        for i in range(len(nmer1_file_lines)):

            if i < 2:
                continue

            else:
                nmer1_atom = nmer1_file_lines[i].split()
                nmer1_info_list.append(nmer1_atom)

        for j in range(len(nmer2_file_lines)):

            if j < 2:
                continue

            else:
                nmer2_atom = nmer2_file_lines[j].split()
                nmer2_info_list.append(nmer2_atom)

        for atom1 in range(len(nmer1_info_list)):

            for atom2 in range(len(nmer2_info_list)):
                x12 = float(nmer1_info_list[atom1][1]) - float(nmer2_info_list[atom2][1])
                y12 = float(nmer1_info_list[atom1][2]) - float(nmer2_info_list[atom2][2])
                z12 = float(nmer1_info_list[atom1][3]) - float(nmer2_info_list[atom2][3])
                r12 = math.sqrt(x12**2.0 + y12**2.0 + z12**2.0)
                atomic_separations.append(r12)

    shortest_atomic_separation = min(atomic_separations)

    return shortest_atomic_separation
# ==================================================================

# ==================================================================
def nmer_merger(nm_a, nm_b, nm_new):
    """Takes the .xyz filenames of two existing N-mers and of a new N-mer to be written.
    Produces the .xyz file of the new N-mer."""

    with open(nm_a, 'r') as nm_a_f, open(nm_b, 'r') as nm_b_f:
        
        new_nm_list = []

        i = 0
        for line_nm_a_f in nm_a_f.readlines():
            i += 1

            if i == 1:
                new_nm_list_line = str(int(line_nm_a_f) + int(nm_b_f.readline())) + "\n"
                new_nm_list.append(new_nm_list_line)

            elif line_nm_a_f.startswith("Monomer"):
                new_nm_list_new_line = "Monomers " + str(nm_a)[2:-4] + "+" + str(nm_b)[2:-4] + "\n"
                new_nm_list.append(new_nm_list_new_line)

            else:
                new_nm_list.append(line_nm_a_f)
        
        new_nm_list.append("\n")

        j = 0
        for line_nm_b_f in nm_b_f.readlines():
            j += 1

            if j < 2:
                continue

            else:
                new_nm_list.append(line_nm_b_f)

    with open(nm_new, "w") as nm_n_f:

        for line_new_nm_list in new_nm_list:
            
            if line_new_nm_list.rstrip():

                nm_n_f.write(line_new_nm_list)

# ==================================================================

# ==================================================================
def nmer_builder(nmer_type, nmer_separation_cutoff):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter."""

    # Function reused for different types of N-mers.

    if nmer_type == "dimers":
        print("Merging monomers with monomers to obtain dimers.\n")
        nm_filename_pattern = "1-*.xyz"
        nm_txt_lbl = "Dimer"
        nm_num_lbl = "2"

    elif nmer_type == "trimers":
        print("Merging dimers with monomers to obtain trimers.\n")
        nm_filename_pattern = "2-*+*.xyz"
        nm_txt_lbl = "Trimer"
        nm_num_lbl = "3"

    elif nmer_type == "tetramers":
        print("Merging trimers with monomers to obtain tetramers.\n")
        nm_filename_pattern = "3-*+*+*.xyz"
        nm_txt_lbl = "Tetramer"
        nm_num_lbl = "4"

    elif nmer_type == "pentamers":
        print("Merging tetramers with monomers to obtain pentamers.\n")
        nm_filename_pattern = "4-*+*+*.xyz"
        nm_txt_lbl = "Pentamer"
        nm_num_lbl = "5"
    
    else:
        print_error("The N-mer type must be defined as 'dimers', 'trimers', 'tetramers', or 'pentamer'.")

    directory = os.getcwd()
    nmer_files = os.listdir(directory)
    num_monomers_ctrl_cell = 12 # TODO: Number of fragments in the central cell for monomer-in-the-cell cutoff.
    counter_monomers  = 0
    counter_nmers     = 0
    counter_new_nmers = 0
    counter_dscrd_spi = 0
    counter_dscrd_exs = 0
    counter_dscrd_sep = 0

    for monomer in nmer_files:
        
        if (re.match('^1-[0-9]+.xyz$', monomer) and counter_monomers < num_monomers_ctrl_cell): # Dimer with at least one monomer in the central cell filter

            for nmer in nmer_files:

                if fnmatch.fnmatch(nmer, nm_filename_pattern):
                    new_nm_filename = nm_num_lbl + "-" + str(monomer)[2:-4] + "+" + str(nmer)[2:]

                    if monomer <= nmer:
                        counter_nmers += 1
                        
                        min_separation = r_min(nmer, monomer)

                        if min_separation <= 0.001: # Superimposed monomers filter. For trimers and tetramers.
                            print( "%s %i (%s) discarded: Separation %3.2f A, superimposed structures" % \
                                   (nm_txt_lbl, counter_nmers, new_nm_filename, min_separation) )
                            counter_dscrd_spi += 1
                            continue

                        elif nmer_separation_cutoff <= min_separation: # Separation cutoff filter.
                            print( "%s %i (%s) discarded: Separation %3.2f A longer than cutoff %3.2f A" \
                                   % (nm_txt_lbl, counter_nmers, new_nm_filename, min_separation, nmer_separation_cutoff) )
                            counter_dscrd_sep += 1

                            # TODO: Nuclear repulsion energy filter.
                            # import psi4
                            # g = psi4.geometry("""He 0 0 0\nHe 0 0 2""")
                            # g.update_geometry()
                            # g.nuclear_repulsion_energy()

                            
                            # TODO: Filter for symetric structures.
 
                        else:
                            nmer_merger(monomer, nmer, new_nm_filename)
                            print( "%s %i (%s) generated: Merged %s and %s" \
                                   % (nm_txt_lbl, counter_nmers, new_nm_filename, monomer, nmer) )
                            counter_new_nmers += 1
                    
                    else: # Double-counting filter
                        counter_nmers += 1
                        print( "%s %i (%s) discarded: %s of %s and %s already exists" \
                               % (nm_txt_lbl, counter_nmers, new_nm_filename, nm_txt_lbl, nmer, monomer) )
                        counter_dscrd_exs += 1

            if (nmer_type == "dimers"):
                counter_monomers += 1		

    print("\nTotal number of monomers: %i" % counter_total_monomers)
    print("\nMaximum possible number of dimers: %i" % counter_total_dimers)
    print("Maximum possible number of trimers: %i" % counter_total_trimers)
    print("Maximum possible number of tetramers: %i" % counter_total_tetramers)
    print("Maximum possible number of pentamers: %i" % counter_total_pentamers)
    print("\nDiscarded %i %s with far-separated atoms." % (counter_dscrd_sep, nm_txt_lbl.lower() + "s"))
    print("Discarded %i %s with double-counted structures." % (counter_dscrd_exs, nm_txt_lbl.lower() + "s"))
    print("Discarded %i %s with superimposing atoms." % (counter_dscrd_spi, nm_txt_lbl.lower() + "s"))
    print("\nGenerated %i %s.\n" % (counter_new_nmers, nm_txt_lbl.lower() + "s"))
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
    
    read_cif_input = 'Benzene-138K.cif'  # CIF input file. # WARNING: Hardcoded for now!
    read_cif_output = 'Benzene-138K.xyz' # XYZ output file. # WARNING: Hardcoded for now!
    read_cif_a = '5'                  # X replicas. # WARNING: Hardcoded for now!
    read_cif_b = '5'                  # Y replicas. # WARNING: Hardcoded for now!
    read_cif_c = '5'                  # Z replicas. # WARNING: Hardcoded for now!
    
    read_cif_arguments = ['', '-i', read_cif_input, '-o', read_cif_output, '-b', read_cif_a, read_cif_b, read_cif_c]
    
    print("The following arguments will be passed to the CIF reader script:\n")
    print("./Read_CIF.py" + " ".join(str(read_cif_argument) for read_cif_argument in read_cif_arguments) + "\n")
    
    Read_CIF.main(read_cif_arguments)
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Take the supercell .xyz file
    # And generate .xyz files with all possible monomers.
    print ("")
    
    # Read the lines in the .xyz file recently generated by Read_CIF.py
    #with open(read_cif_output) as supercell_xyz_file:
    #    sc_xyz_f_lines = supercell_xyz_file.readlines()
    
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

    fragmented_supecell_file = "bztest.p4"      # NOTE: for test only!
    #fragmented_supecell_file = "bzfrag.p4"     # TODO: Name of the fragmented super cell file is temporarily hardcoded.
    
    num_frags = 664           # TODO: Number of fragments temporarily hardcoded!
    num_atoms_frag = 12       # TODO: Number of atoms per fragment hardcoded! What if there are two types of molecules?
    fragment_separator = "--" # Fragment separator string.

    psifragments2xyzs(fragmented_supecell_file, num_frags, num_atoms_frag, fragment_separator)

    # Discard fragments that are not a complete molecule.
    
    print("Detecting fragments with incomplete molecules.\n")

    fragment_filter(num_atoms_frag)

    # Organize fragments according to their distance to the center of
    # the supercell, produce molecule files.
    
    print("Organizing molecules according to their separation to the center of the supercell.\n")

    fragments2molecules()
    
    # Create monomer files from fragment files.

    print("Creating monomers from molecule files.\n")

    molecules2monomers()
    
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all monomers and generate dimers with all other
    # monomers.
    #
    # Filter dimers that are too distant apart.
    
    r_cut_dimer = 3.0
    nmer_builder("dimers", r_cut_dimer)

    # Filter out and keep count of all non-unique dimers, using the
    # nuclear repulsion energy criteria.
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all dimers and generate trimers with all other
    # monomers.
    #
    # Filter trimers that are too distant apart.

    r_cut_trimer = 3.0
    nmer_builder("trimers", r_cut_trimer)

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
