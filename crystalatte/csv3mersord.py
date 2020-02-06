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

import csv
import operator
import os

def reorder_file(fname):

    with open(fname, "r") as csvf:
    
        # Read the CSV file and get the data
        csvdata = csv.reader(csvf, delimiter=",")
        
        # Skip first line
        next(csvdata) 
    
        # Parse each row into the proper python type.
        parsed = [[row[0],
                   float(row[1]), 
                   int(row[2]), 
                   float(row[3]), 
                   float(row[4]), 
                   float(row[5]), 
                   float(row[6]), 
                   float(row[7]), 
                   float(row[8])]
                  for row in csvdata]
    
        # Sort by contribution and then by each separation.
        sortset1 = sorted(parsed, key=operator.itemgetter(3), reverse=True)
        sortset2 = sorted(sortset1, key=operator.itemgetter(6))
        sortset3 = sorted(sortset2, key=operator.itemgetter(7))
        sortsetf = sorted(sortset3, key=operator.itemgetter(8))
    
        # The crystal lattice energy will be recomputed.
        cle = 0.0
    
        # Generate new file name.
        try:
            newcsv = os.path.splitext(fname)[0] + "-3mC.csv"
        except TypeError:
            newcsv = "Results.csv"
        except NameError:
            newcsv = "Results.csv"
    
        # Write the new file, recomputing the accumulated crystal lattice energy.
        with open(newcsv, 'w') as csvf:
    
            csv_header = "N-mer Name,"\
                    + "Non-Additive MB Energy (kJ/mol),"\
                    + "Num. Rep. (#), N-mer Contribution (kJ/mol),"\
                    + "Partial Crys. Lattice Ener. (kJ/mol),"\
                    + "Calculation Priority (Arb. Units),"\
                    + "Minimum Monomer Separations (A)"
    
            csvf.write(csv_header + "\n")
        
            for row in sortsetf:
    
                # Recompute the crystal lattice energy.
                cle += row[3]
    
                # Format the string for each line in the new CSV file.
                nmer_csv = "{:},{:.8f},{:},{:.8f},{:.8f},{:.6e},{:.3f},{:.3f},{:.3f}".format(
                        row[0], row[1], row[2], row[3], cle, row[5], row[6], row[7], row[8])
    
                # Write the line on the new CSV file.
                csvf.write(nmer_csv + "\n")
    
    return newcsv

def orderer_main():

    d = os.getcwd()

    for f in os.listdir(d):

        if f.endswith(".csv"):
            reorder_file(f)

        else:
            continue

    return

if __name__ == "__main__":
    orderer_main()
