#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice
# energies.
#
# Copyright (c) 2017-2018
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


# ======================================================================
def success_check(fname, verbose=0):
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
def get_nmer_data(fname, verbose=0):
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

        n_body_energy = None

        for line in outf:

            # Increment the line counter by 1.
            i += 1

            # Find the name of the N-mer, and use it to create
            # the key for its entry in the nmers dictionary.
            if "Psithon input for N-mer:" in line:
                splt = line[:-1].split(":")
                keynmer = splt[-1].strip()
            
            # Find the number of replicas of the N-mer.
            if "Number of replicas:" in line:
                splt = line[:-1].split(":")
                replicas = int(splt[-1].strip())

            # Find the number of replicas of the N-mer.
            if "Priority index for input:" in line:
                splt = line[:-1].split(":")
                priority_min = float(splt[-1].strip())

            # Find the list of minimum monomer separations of the N-mer.
            if "Minimum monomer separations:" in line:
                splt = line[:-1].split(":")
                msps = splt[-1].strip()
                min_monomer_separations = msps.split()

            # Find the list of COM monomer separations of the N-mer.
            if "Minimum COM separations:" in line:
                splt = line[:-1].split(":")
                csps = splt[-1].strip()
                com_monomer_separations = csps.split()
            
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
                # number is in Kcal/mol and is now converted to KJ/mol
                elif not nbend and i > nbini:
                    lastl = line.strip().split()
                    number_of_monomers = int(lastl[0])
                    n_body_energy = float(lastl[-1]) * 4.184 # Same value as in qcdb.psi_cal2J
                    #print(line[:-1]) #debug
            
    return keynmer, number_of_monomers, replicas, priority_min, min_monomer_separations, com_monomer_separations, n_body_energy
    
# ======================================================================

# ======================================================================
def print_results(results, crystal_lattice_energy, verbose=0):
    """Prints a summary of the energy results at the end of the
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
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------------------------------------------------------------")
        print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (KJ/mol) |  (#) |     (KJ/mol) |      (KJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------------------------------------------------------------")
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------------------------------------------------------------")
        print("\nCrystal Lattice Energy (Eh)       = {:5.8f}".format(crystal_lattice_energy * 2625.500)) # Same value as in psi_hartree2kJmol
        print("Crystal Lattice Energy (KJ/mol)   = {:9.8f}".format(crystal_lattice_energy))
        print("Crystal Lattice Energy (Kcal/mol) = {:9.8f}\n".format(crystal_lattice_energy / 4.184)) # Same value as in psi_cal2J
# ======================================================================


# ======================================================================
def print_end_msg(start, verbose=0):
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
def main(verbose=0):

    # Start counting execution time.
    start = time.time()

    d = os.getcwd()
    
    # A dictionary for storing N-mers will be created.
    nmers = {}

    crystal_lattice_energy = 0.0
    results = []
    csv_lines = []
    
    csv_header = "N-mer Name,"\
            + "Non-Additive MB Energy (KJ/mol),"\
            + "Num. Rep. (#), N-mer Contribution (KJ/mol),"\
            + "Partial Crys. Lattice Ener. (KJ/mol),"\
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
            success = success_check(f)
            #print(success) #debug
    
            if success:
                # If the output was ran successufully, get the N-mer name
                # to create a key for its soon-to-be created dictionary
                # and get the data to populate such N-mer dictionary.
                keynmer, number_of_monomers, replicas, priority_min, min_monomer_separations, com_monomer_separations, n_body_energy = get_nmer_data(f)

                # Create a new dictionary for the current N-mer.
                nmers[keynmer]= {}

                # Populate the current N-mer dictionary with the data
                # from the Psi4 output file.
                nmers[keynmer]["replicas"] = replicas
                nmers[keynmer]["priority_min"] = priority_min
                nmers[keynmer]["min_monomer_separations"] = min_monomer_separations
                nmers[keynmer]["com_monomer_separations"] = com_monomer_separations
                nmers[keynmer]["nambe"] = n_body_energy 
                nmers[keynmer]["contrib"] = n_body_energy * replicas / number_of_monomers
                
                crystal_lattice_energy += nmers[keynmer]["contrib"]

                #print(n_body_energy) #debug
                #print(replicas) #debug
                #print(number_of_monomers) #debug
                #print(nmers[keynmer]["contrib"]) #debug
                #print(crystal_lattice_energy) #debug
                
                # Generate a string with an ordered list of minimum separations
                # between atoms belonging to different monomers.
                rminseps = ""
                
                nmer_min_monomer_separations = nmers[keynmer]["min_monomer_separations"] 
                nmer_min_monomer_separations.sort()
                
                for r in nmer_min_monomer_separations:
                    rminseps += "{:6.3f} ".format(float(r))

                nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                        keynmer,
                        nmers[keynmer]["nambe"],
                        nmers[keynmer]["replicas"],
                        nmers[keynmer]["contrib"],
                        crystal_lattice_energy,
                        nmers[keynmer]["priority_min"],
                        rminseps)
                
                #print(nmer_result) #debug
                results.append(nmer_result)

                nmer_csv = "{:26} , {:>12.8f} , {:>4} , {:>12.8f} , {:>13.8f} , {:12.6e} , {}".format(
                        keynmer,
                        nmers[keynmer]["nambe"],
                        nmers[keynmer]["replicas"],
                        nmers[keynmer]["contrib"],
                        crystal_lattice_energy,
                        nmers[keynmer]["priority_min"],
                        rminseps)
                
                #print(nmer_csv) #debug
                csv_lines.append(nmer_csv) #debug

        else:
            continue
    
    print_results(results, crystal_lattice_energy, verbose)
    print_end_msg(start, verbose)

    #TODO: This should be the CIF file name, when that is inclued in 
    #      the Psi4 output
    csvname = "Results.csv"

    with open(csvname, 'w') as csvf:
        for line in csv_lines:
            #print(line) #debug
            csvf.write(line + "\n")
    
    #print("\nThe N-mers Dictionary looks like this:\n") #debug
    #print(nmers) #debug

if __name__ == "__main__":
    main(1)
