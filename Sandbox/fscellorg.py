#!/usr/bin/env python

import math
import re
import os

# Read the coordinates of all fragments and find the origin.

def readmatchf(pattern, file):
    """Read lines contained in a given file if its name matches a pattern."""
    
    if re.match(pattern, file):  # Match filenames with pattern
        
        with open(file, 'r') as f:
            contents = f.readlines()

    return contents

def main():
    """Compute the center of coordinates of a series of .xyz files with a matching pattern in the filename."""
    
    directory = os.getcwd()
    files = os.listdir(directory)

    x = []       # Values of x coordinates.
    y = []       # Values of y coordinates.
    z = []       # Values of z coordinates.
    data = {}    # File names and text inside.
    txt = []     # Lines with coordinates.
    skipped = 0  # Counts number of lines with no atom coordinates.


    for file in files:
        
        if re.match('^f[0-9]+.xyz$', file):  # Match filenames f???.xyz
            
            with open(file, 'r') as f:
                lines = f.readlines()

                i = 0

                for i in range(len(lines)):

                    if i < 2:
                        skipped += 1
                    
                    else:
                        txt.append(lines[i])
                        info = lines[i].split()
                        x.append(float(info[1]))
                        y.append(float(info[2]))
                        z.append(float(info[3]))

    xavg = (max(x) - min(x))/2.0
    yavg = (max(y) - min(y))/2.0
    zavg = (max(z) - min(z))/2.0

    r = []

    for j in range(len(x)):
        r.append(math.sqrt((x[j] - xavg)**2.0 + (y[j] - yavg)**2.0 + (z[j] - zavg)**2.0))

    print(min(r))
    print(r.index(min(r))) # Index of the atom closest to the origin.
    print(txt[r.index(min(r))]) # Line with the atom closest to the origin.

    return

if __name__ == "__main__":
    main()
