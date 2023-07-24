#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice
# energies.
#
# Copyright (c) 2017-2020
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
import os, sys
import time
import shutil
import re
import numpy as np

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
        txt = outf.read()
        match  = re.search("Buy a developer a beer!", txt)
        if match:
            beer = True

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

        lines = outf.readlines()
        txt = ''.join(lines[:250])

        if (match := re.search(r"^.*# Psithon input for N-mer:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            key = splt[-1].strip()

        if (match := re.search(r"^.*# Generated from:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            xyz_path = splt[-1].strip()
            tree = xyz_path.split("/")
            sc = tree[-1].strip()
            # Remove .xyz extension in case it has it.
            if sc.endswith(".xyz"):
                sc_xyz = sc[:-4]
            else:
                sc_xyz = sc
                   
        # Find the number of replicas of the N-mer.
        if (match := re.search(r"^.*# Number of replicas:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            replicas = int(splt[-1].strip())

        # Find the list of minimum monomer separations of the N-mer.
        if (match := re.search(r"^.*# Minimum monomer separations:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
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
        if (match := re.search(r"^.*# Minimum COM separations:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            csps = splt[-1].strip()
            com_monomer_separations = csps.split()

        # NOTE: Strings deprecated in nmer2psithon() function.
        # These chunks are kept here for backward compatibility.
        if (match := re.search(r"^.*# COM priority index for input:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            priority_com = float(splt[-1].strip())

        if (match := re.search(r"^.*# Priority index for input:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            priority_min = float(splt[-1].strip())

        # Get the cutoff-based priority of the N-mer.
        if (match := re.search(r"^.*# Cutoff priority:.*$",txt, re.MULTILINE)):
            splt = match.group().split(":")
            priority_cutoff = float(splt[-1].strip())

        # Get the atomic-separation-based priority of the N-mer.
        if (match := re.search(r"^.*# Separation priority:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            priority_min = float(splt[-1].strip())

        # Get the COM-separation-based priority of the N-mer.
        if (match := re.search(r"^.*# COM priority:.*$", txt, re.MULTILINE)):
            splt = match.group().split(":")
            priority_com = float(splt[-1].strip())

        # Find the nuclear repulsion energy of the N-mer.
        if (match := re.search(r"^.*# Nuclear repulsion energy:.*$", txt, re.MULTILINE)):

            # Remove units if present.
            if "a.u." in match.group().lower():
                text = match.group()[:-5]

            else:
                text = match.group()

            splt = text.split(":")
            nre = float(splt[-1].strip())

        txt = ''.join(lines[-250:])

        if (match := re.search(r"n-Body.*?\n{2,}?", txt, re.S)):
            lastl = match.group().strip().split('\n')[-1]
            lastl = lastl.split()
            if lastl[0] == "FULL/RTN":
                lastl = lastl[1:]
            number_of_monomers = int(lastl[0])
            n_body_energy = float(lastl[-1]) * 4.184 # Same value as in qcdb.psi_cal2J

        # If it's a SAPT computation, we grab the 2-body energy in
        # a different output line
        if (match := re.search(r"^.*Total SAPT.*$", txt, re.MULTILINE)):
            n_body_energy = float(match.group().split()[-2])
            number_of_monomers = 2

    # Just in case the COM priority was not defined.
    try:
        priority_com

    except UnboundLocalError:
        priority_com = priority_min

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
def psz_print_results(results, crystal_lattice_energy, com_mode, verbose=0):
    """Prints a summary of the energy results at the end of the
    execution.
    
    Arguments:
    <list> results
        List of results
    <float> crystal_lattice_energy
        Total crystal lattice energy (kJ/mol)
    <bool> com_mode
        Print center of mass distances instead of interatomic distances
    <int> verbose
        Adjusts the level of detail of the printouts.
    """
    try:
        term_size = shutil.get_terminal_size().columns
    except AttributeError:
        term_size = 0

    if verbose >= 1:
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(term_size - 112)))
        if not com_mode:
            print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        else:
            print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | COM Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (kJ/mol) |  (#) |     (kJ/mol) |      (kJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(term_size - 112)))
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}\n".format("-"*(term_size - 112)))
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
def psz_main(verbose=0, **kwargs):
    """Main routine to process psithon output files, print a summary, and dump
    corresponding .dat and .csv files.

    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    **kwargs
        Keyword arguments as specified in psz_process_args_and_run()
    """

    # Start counting execution time.
    start = time.time()

    if kwargs['src_directory'] == "":
        d = os.getcwd()
    else:
        d = kwargs['src_directory']
    
    com_mode = kwargs['com_mode']

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
            + "Calculation Priority (Arb. Units)," \
            + "Avg COM Separation (A)," \
            + "Avg Monomer Separation (A)," \

    if not com_mode:
        csv_header = csv_header + "Minimum Monomer Separations (A)"
    else:
        csv_header = csv_header + "COM Monomer Separations (A)"
        

    csv_lines.append(csv_header)
    #print(csv_header) #degub

    for f in os.listdir(d):
    
        # Find output files using the default extension.
        #NOTE: This is prone to error if the extension was changed.
        if f.endswith(".out"):
            
            #print(os.path.join(d, f)) #debug

            filename = os.path.join(d, f)

            # Check if Psi4 exited successfully.
            # Why bothering analyzing a file otherwise?
            success = psz_success_check(filename)
            #print(success) #debug
    
            if success:
                # If the output was ran successufully, get the N-mer name
                # to create a key for its soon-to-be created dictionary
                # and get the data to populate such N-mer dictionary.
                sc_xyz, key, number_of_monomers, replicas, p_cutoff, p_min, p_com, min_mon_seps, com_mon_seps, nre, n_body_energy = psz_get_nmer_data(filename)

                # Create a new dictionary for the current N-mer.
                nmers[key]= {}

                # Populate the current N-mer dictionary with the data
                # from the Psi4 output file.
                nmers[key]["replicas"] = replicas
                nmers[key]["priority_cutoff"] = p_cutoff
                nmers[key]["priority_min"] = p_min
                nmers[key]["priority_com"] = p_com
                nmers[key]["min_monomer_separations"] = min_mon_seps
                nmers[key]["average_min_monomer_separation"] = np.average(np.array(min_mon_seps, dtype=np.float64))
                nmers[key]["com_monomer_separations"] = com_mon_seps
                nmers[key]["average_com_monomer_separation"] = np.average(np.array(com_mon_seps, dtype=np.float64))
                nmers[key]["nre"] = nre
                nmers[key]["nambe"] = n_body_energy 
                nmers[key]["contrib"] = n_body_energy * replicas / number_of_monomers
                
                crystal_lattice_energy += nmers[key]["contrib"]

        else:
            continue

    # Get the keys of the N-mers dictionary, and put them on a list.
    nmer_keys = list(nmers.keys())

    # Sort the list in decreasing priority order.
    if kwargs['sort_by_avg_com_dist']:
        nmer_keys.sort(key = lambda x: nmers[x]['average_com_monomer_separation'])
    elif kwargs['sort_by_nmer_cutoff']:
        nmer_keys.sort(key = lambda x: float(nmers[x]['min_monomer_separations'][-1]))
    else:
        nmer_keys.sort(key = lambda x: -nmers[x]['priority_cutoff'])

    # The next line was replaced to trigger the calculations in order.
    for keynmer in nmer_keys:
        #nmer = nmers[keynmer]

        partial_crystal_lattice_energy += nmers[keynmer]["contrib"]

        # Generate a string with an ordered list of minimum separations
        # between atoms belonging to different monomers.
        rminseps    = ""
        rminsepscsv = ""
        
        # Same as above but for COM separations
        comseps    = ""
        comsepscsv = "" 
        
        nmer_min_monomer_separations = nmers[keynmer]["min_monomer_separations"] 
        
        for r in nmer_min_monomer_separations:
            rminseps    += "{:6.3f} ".format(float(r))
            rminsepscsv += "{:.3f},".format(float(r))

        for r in nmers[keynmer]["com_monomer_separations"]:
            comseps    += "{:6.3f} ".format(float(r))
            comsepscsv += "{:.3f},".format(float(r))

        nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                keynmer,
                nmers[keynmer]["nambe"],
                nmers[keynmer]["replicas"],
                nmers[keynmer]["contrib"],
                partial_crystal_lattice_energy,
                nmers[keynmer]["priority_cutoff"],
                rminseps if not com_mode else comseps)
        
        results.append(nmer_result)

        nmer_csv = "{:},{:.8f},{:},{:.8f},{:.8f},{:.6e},{:.4f},{:.4f},{}".format(
                keynmer,
                nmers[keynmer]["nambe"],
                nmers[keynmer]["replicas"],
                nmers[keynmer]["contrib"],
                partial_crystal_lattice_energy,
                nmers[keynmer]["priority_min"],
                nmers[keynmer]["average_com_monomer_separation"],
                nmers[keynmer]["average_min_monomer_separation"],
                rminsepscsv[:-1] if not com_mode else comsepscsv[:-1])
        
        csv_lines.append(nmer_csv)

    psz_print_header(verbose)
    psz_print_results(results, crystal_lattice_energy, com_mode, verbose)
    psz_print_end_msg(start, verbose)

    try:
        csvname = sc_xyz + ".csv"

    except TypeError:
        csvname = "Results.csv"

    except NameError:
        csvname = "Results.csv"

    with open(csvname, 'w') as csvf:

        for line in csv_lines:
            csvf.write(line + "\n")

    return results, crystal_lattice_energy
    

def psz_process_args_and_run():
    """Process command line arguments and run the script.

    psithonizer.py [--com_mode] [--sort_by_avg_com_dist] [src_directory]

    --com_mode
        print center of mass separations rather than minimum monomer separations
        (default False)
    --sort_by_nmer_cutoff
        sort by the nmer cutoff (largest of the intermonomer separations)
    --sort_by_avg_com_dist
        sort the output lines by the average COM distance
        (default False)
    src_directory
        directory of the output files to analyze (defaults to current working directory)
    """

    kwargs = {'com_mode': False, 'sort_by_avg_com_dist': False, 'sort_by_nmer_cutoff': False, 'src_directory': ""}

    arglist = sys.argv
    arglist.pop(0) # delete first argument, that is the name of this script

    for arg in arglist:
        if (arg == "--com_mode"):
            kwargs['com_mode'] = True
        elif (arg == "--sort_by_avg_com_dist"):
            kwargs['sort_by_avg_com_dist'] = True
        elif (arg == "--sort_by_nmer_cutoff"):
            kwargs['sort_by_nmer_cutoff'] = True
        else: # assume we are giving a directory to some output files
            kwargs['src_directory'] = arg

    if (kwargs['sort_by_avg_com_dist'] and kwargs['sort_by_nmer_cutoff']):
        print("Can only sort by one thing at a time!")
        exit()

    psz_main(1, **kwargs)


if __name__ == "__main__":

    psz_process_args_and_run()

