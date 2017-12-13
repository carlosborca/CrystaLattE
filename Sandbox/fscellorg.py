#!/usr/bin/env python

import math
import re
import os

# Read the coordinates of all fragments and find the origin.

def readmatchf(pattern, file):
    """If the name of a file matches a pattern read the lines contained in it, and return the contents as a list of lines."""

    contents = []

    if re.match(pattern, file):  # Match filenames with pattern
        
        with open(file, 'r') as f:
            contents = f.readlines()

    return contents


def scellcntr():
    """Compute the center of coordinates of a series of .xyz files with a matching pattern in the filename.
    Returns a list with fragment names, x, y, and z coordinates in lists of floats and the coordinates of the center."""
    
    frg = [] # Fragments file names.
    x = [] # Values of x coordinates.
    y = [] # Values of y coordinates.
    z = [] # Values of z coordinates.

    directory = os.getcwd()
    files = os.listdir(directory)

    for file in files:
        lines = readmatchf('^f[0-9]+.xyz$', file)
        
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

    return frg, x, y, z, cntr_x, cntr_y, cntr_z

def main():

    frg, x, y, z, cntr_x, cntr_y, cntr_z = scellcntr()

    r = []

    for j in range(len(x)):
        r.append(math.sqrt((x[j] - cntr_x)**2.0 + (y[j] - cntr_y)**2.0 + (z[j] - cntr_z)**2.0))

    print(min(r))
    print(r.index(min(r))) # Index of the atom closest to the origin.
    print(frg[r.index(min(r))]) # Name of the file with the atom closests to the origin.

    return

if __name__ == "__main__":
    main()
