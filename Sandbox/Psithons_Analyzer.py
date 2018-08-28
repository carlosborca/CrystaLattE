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
import subprocess

# A dictionary for storing N-mers will be created.
nmers = {}

d = os.getcwd()

for f in os.listdir(d):

    # Find output files using the default extension.
    # WARNING: This is prone to error if the extension was changed.
    if f.endswith(".out"):
        
        print(os.path.join(d, f))

        # Check if Psi4 exited successfully.
        # Why bothering analyzing a file otherwise?
        beer = False
        with open(f, 'r') as outf:

            for line in outf:
                if "Psi4 exiting successfully. Buy a developer a beer!" in line:
                    beer = True
                    break
                else:
                    continue
        
        # If the output was ran successufully, get the data.
        if beer:
            with open(f, 'r') as outf:
                print(f)
                for line in outf:
                    if "Psithon input for N-mer:" in line:
                        print(line)
            
            print("Success!")
            
        else:
            print("Failure!")

    else:
        continue
