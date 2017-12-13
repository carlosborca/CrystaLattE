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

def posvecmag():
    """Using the scellcntr() function it computes the magnitude of the position vector of the atoms read by the scellcntr() function.
    Returns the fragment filename and the magnitude of the position vector."""

    frg, x, y, z, cntr_x, cntr_y, cntr_z = scellcntr()

    r = []

    for j in range(len(x)):
        r.append(math.sqrt((x[j] - cntr_x)**2.0 + (y[j] - cntr_y)**2.0 + (z[j] - cntr_z)**2.0))

    return frg, r

def proximity():
    """Using position vector magnitude, computed in posvecmag(). This function returns a list of fragments filenames
    with fragments containing atoms closest to the center of the supercell first."""

    frg, r = posvecmag()
    
    frgprox = []

    for pair in sorted(enumerate(r), key=lambda item:item[1]):
        if frg[pair[0]] in frgprox:  
            continue
        else:
            frgprox.append(frg[pair[0]])

    return frgprox

def main():
    """Main program."""
    
    directory = os.getcwd()
    proxlist = proximity()
    molecule = 0

    for file in proxlist:
        molecule += 1
        molfname = "m" + str(molecule).zfill(len(proxlist[-1]) - 5) + ".xyz"

        with open(file, 'r') as f, open(molfname, 'w') as m:
            
            for line in file:

                if line.startswith("Fragment"):
                    newline = "Molecule " + str(molecule) + " (" + line[:-1] + ")"
                    m.write(newline)

                else:
                    m.write(line)
        
        os.remove(os.path.join(directory,file))

    return

if __name__ == "__main__":
    main()
