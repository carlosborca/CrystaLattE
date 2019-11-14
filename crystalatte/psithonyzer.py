#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice
# energies.
#
# Copyright (c) 2017-2019
# Carlos H. Borca
# C. David Sherrill
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
import os
import time
import shutil

# ======================================================================
def psz_success_check(fname, verbose=0):
    """Takes a Psi4 output file and checks if it contains the string
    'coffee' to decide if the calculation was successful or not.

    Arguments:
    <str> fname
        Name of the Psi4 output file to check.
    <int> verbose
        Adjusts the level of detail of the printouts.
    
    Returns:
    <bool> beer
        True if the calculation was successful, false otherwise.
    """

    beer = False
    with open(fname, 'r') as outf:
        for line in outf:
            if "Psi4 exiting successfully. Buy a developer a beer!" in line:
                beer = True
                break
            
            else:
                continue
    return beer
# ======================================================================


# ======================================================================
def psz_get_nmer_data(fname, verbose=0):
    """.
    
    Arguments:
    <str> fname
        Name of the Psi4 output file to check.
    <int> verbose
        Adjusts the level of detail of the printouts.
    
    Returns:
    <str> key
        The key to create a new N-mer dictionary in the dictionary
        containing all N-mers.
    .
    .
    .

    """

    with open(fname, 'r') as outf:
        
        # Start a line counter
        i = 0
        
        # Index to indicate where the N-body decomposition
        # information starts and ends in the Psi4 output file.
        nbini = None
        nbend = None

        # Initialize the value of required variables.
        sc_xyz = None
        nre = None
        n_body_energy = None

        for line in outf:

            # Increment the line counter by 1.
            i += 1

            # Find the name of the N-mer, and use it to create
            # the key for its entry in the nmers dictionary.
            if "Psithon input for N-mer:".lower() in line.lower():
                splt = line[:-1].split(":")
                key = splt[-1].strip()

            # NOTE: This only works on Linux and MacOS.
            # Find the name of the XYZ from which this file generated.
            if "Generated from:".lower() in line.lower():
                splt = line[:-1].split(":")
                xyz_path = splt[-1].strip()
                tree = xyz_path.split("/")
                sc = tree[-1].strip()

                # Remove .xyz extension in case it has it.
                if sc.endswith(".xyz"):
                    sc_xyz = sc[:-4]

                else:
                    sc_xyz = sc

            # Find the number of replicas of the N-mer.
            if "Number of replicas:".lower() in line.lower():
                splt = line[:-1].split(":")
                replicas = int(splt[-1].strip())

            # Find the list of minimum monomer separations of the N-mer.
            if "# Minimum monomer separations:".lower() in line.lower():
                splt = line[:-1].split(":")
                msps = splt[-1].strip()
                min_monomer_separations = msps.split()

                # NOTE: Because the cutoff priority was not originally
                # implemented in CrystaLattE, if it is not found in the
                # psithon output, it will be approximated from the list
                # of minimum monomer separations.
                max_sep = float(max(min_monomer_separations))

                if len(min_monomer_separations) == 1:
                    num_mon = 2.0

                if len(min_monomer_separations) == 3:
                    num_mon = 3.0

                if len(min_monomer_separations) == 6:
                    num_mon = 4.0

                if len(min_monomer_separations) == 10:
                    num_mon = 5.0

                main_contrib = 1.0/(max_sep**(num_mon**2.0))

                idx = 0
                add = 0.0
                for rmin in min_monomer_separations:

                    if idx != 0:
                        one_over_rmin3 = 1.0e-10/(float(rmin)**(num_mon**2.0))
                        add += one_over_rmin3

                    idx += 1

                priority_cutoff = main_contrib + add

            # Find the list of COM monomer separations of the N-mer.
            if "# Minimum COM separations:".lower() in line.lower():
                splt = line[:-1].split(":")
                csps = splt[-1].strip()
                com_monomer_separations = csps.split()

            # NOTE: Strings deprecated in nmer2psithon() function.
            # These chunks are kept here for backward compatibility.
            if "# COM Priority index for input:".lower() in line.lower():
                splt = line[:-1].split(":")
                priority_com = float(splt[-1].strip())

            if "# Priority index for input:".lower() in line.lower():
                splt = line[:-1].split(":")
                priority_min = float(splt[-1].strip())

            # Get the cutoff-based priority of the N-mer.
            if "# Cutoff priority".lower() in line.lower():
                splt = line[:-1].split(":")
                priority_cutoff = float(splt[-1].strip())

            # Get the atomic-separation-based priority of the N-mer.
            if "# Separation priority".lower() in line.lower():
                splt = line[:-1].split(":")
                priority_min = float(splt[-1].strip())

            # Get the COM-separation-based priority of the N-mer.
            if "# COM priority".lower() in line.lower():
                splt = line[:-1].split(":")
                priority_com = float(splt[-1].strip())

            # Find the nuclear repulsion energy of the N-mer.
            if "# Nuclear repulsion energy:".lower() in line.lower():
                splt = line[:-1].split(":")
                nre = float(splt[-1].strip())

            # Find where the start of the N-Body decomposition 
            # information is.
            if "n-Body" in line:
                nbini = i

            # Find where the end of the N-Body decomposition 
            # information is.
            if nbini:

                # If the line is blank, record the line number as the
                # end of the N-Body decomposition block, unless it has
                # been recorded alredy.
                if not line.strip():
                    if nbend:
                        continue
                    else:
                        nbend = i

                # Take the N-Body energy and the number of monomers in
                # the N-mer from the last line of the N-Body block. This
                # number is in kcal/mol and is now converted to kJ/mol
                elif not nbend and i > nbini:
                    lastl = line.strip().split()
                    number_of_monomers = int(lastl[0])
                    n_body_energy = float(lastl[-1]) * 4.184 # Same value as in qcdb.psi_cal2J

    return sc_xyz, key, number_of_monomers, replicas, priority_cutoff, priority_min, priority_com, min_monomer_separations, com_monomer_separations, nre, n_body_energy

# ======================================================================

# ======================================================================
def psz_print_header(verbose=0):
    """Prints a the header of the program when starting the
    execution.
    
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """
    
    if verbose >= 1:
        print("") 
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("                              CrystaLattE                              \n")
        print("  The tool for the automated calculation of crystal lattice energies.  \n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("CrystaLattE is being executed in Psithon analysis mode.\n")

# ======================================================================

# ======================================================================
def psz_print_results(results, crystal_lattice_energy, verbose=0):
    """Prints a summary of the energy results at the end of the
    execution.
    
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """

    if verbose >= 1:
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (kJ/mol) |  (#) |     (kJ/mol) |      (kJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}\n".format("-"*(shutil.get_terminal_size().columns - 112)))
        #print("Crystal Lattice Energy (Eh)       = {:5.8f}".format(crystal_lattice_energy / 2625.500)) # Same value as in psi_hartree2kJmol
        print("Crystal Lattice Energy (kJ/mol)   = {:9.8f}".format(crystal_lattice_energy))
        print("Crystal Lattice Energy (kcal/mol) = {:9.8f}\n".format(crystal_lattice_energy / 4.184)) # Same value as in psi_cal2J
# ======================================================================


# ======================================================================
def psz_print_end_msg(start, verbose=0):
    """.
    """

    if verbose >= 1:
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("Execution terminated succesfully.")
        print("Total elapsed wall-clock time: {:.2f} s\n".format(time.time() - start))
        print("Thank you for using CrystaLatte.\n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
# ======================================================================


# ======================================================================
def psz_main(verbose=0):

    # Start counting execution time.
    start = time.time()

    d = os.getcwd()
    
    # A dictionary for storing N-mers will be created.
    nmers = {}

    crystal_lattice_energy = 0.0
    partial_crystal_lattice_energy = 0.0

    results = []
    csv_lines = []
    
    csv_header = "N-mer Name,"\
            + "Non-Additive MB Energy (kJ/mol),"\
            + "Num. Rep. (#), N-mer Contribution (kJ/mol),"\
            + "Partial Crys. Lattice Ener. (kJ/mol),"\
            + "Calculation Priority (Arb. Units),"\
            + "Minimum Monomer Separations (A)"

    csv_lines.append(csv_header)
    #print(csv_header) #degub

    for f in os.listdir(d):
    
        # Find output files using the default extension.
        #NOTE: This is prone to error if the extension was changed.
        if f.endswith(".out"):
            
            #print(os.path.join(d, f)) #debug

            # Check if Psi4 exited successfully.
            # Why bothering analyzing a file otherwise?
            success = psz_success_check(f)
            #print(success) #debug
    
            if success:
                # If the output was ran successufully, get the N-mer name
                # to create a key for its soon-to-be created dictionary
                # and get the data to populate such N-mer dictionary.
                sc_xyz, key, number_of_monomers, replicas, p_cutoff, p_min, p_com, min_mon_seps, com_mon_seps, nre, n_body_energy = psz_get_nmer_data(f)

                # Create a new dictionary for the current N-mer.
                nmers[key]= {}

                # Populate the current N-mer dictionary with the data
                # from the Psi4 output file.
                nmers[key]["replicas"] = replicas
                nmers[key]["priority_cutoff"] = p_cutoff
                nmers[key]["priority_min"] = p_min
                nmers[key]["priority_com"] = p_com
                nmers[key]["min_monomer_separations"] = min_mon_seps
                nmers[key]["com_monomer_separations"] = com_mon_seps
                nmers[key]["nre"] = nre
                nmers[key]["nambe"] = n_body_energy 
                nmers[key]["contrib"] = n_body_energy * replicas / number_of_monomers
                
                crystal_lattice_energy += nmers[key]["contrib"]

        else:
            continue

    # Get the keys of the N-mers dictionary, and put them on a list.
    nmer_keys = list(nmers.keys())

    # Sort the list in decreasing priority order.
    nmer_keys.sort(key = lambda x: -nmers[x]['priority_cutoff'])

    # The next line was replaced to trigger the calculations in order.
    for keynmer in nmer_keys:
        #nmer = nmers[keynmer]

        partial_crystal_lattice_energy += nmers[keynmer]["contrib"]

        # Generate a string with an ordered list of minimum separations
        # between atoms belonging to different monomers.
        rminseps    = ""
        rminsepscsv = ""
        
        nmer_min_monomer_separations = nmers[keynmer]["min_monomer_separations"] 
        #nmer_min_monomer_separations.sort() #debug
        
        for r in nmer_min_monomer_separations:
            rminseps    += "{:6.3f} ".format(float(r))
            rminsepscsv += "{:6.3f},".format(float(r))

        nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                keynmer,
                nmers[keynmer]["nambe"],
                nmers[keynmer]["replicas"],
                nmers[keynmer]["contrib"],
                partial_crystal_lattice_energy,
                nmers[keynmer]["priority_cutoff"],
                rminseps)
        
        results.append(nmer_result)

        nmer_csv = "{:26} , {:>12.8f} , {:>4} , {:>12.8f} , {:>13.8f} , {:12.6e} , {}".format(
                keynmer,
                nmers[keynmer]["nambe"],
                nmers[keynmer]["replicas"],
                nmers[keynmer]["contrib"],
                partial_crystal_lattice_energy,
                nmers[keynmer]["priority_min"],
                rminsepscsv)
        
        csv_lines.append(nmer_csv) #debug

    psz_print_header(verbose)
    psz_print_results(results, crystal_lattice_energy, verbose)
    psz_print_end_msg(start, verbose)

    try:
        csvname = sc_xyz + ".csv"

    except NameError:
        csvname = "Results.csv"

    with open(csvname, 'w') as csvf:
        for line in csv_lines:
            csvf.write(line + "\n")
    
if __name__ == "__main__":
    psz_main(1)
