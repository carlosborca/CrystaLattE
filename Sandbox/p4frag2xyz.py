#!/usr/bin/env python

p4frag = "bzfrag.p4" # WARNING: Name of the fragmented super cell file is temporarily hardcoded. Would be provided by CrystaLattE.

with open(p4frag) as cellfrags:
    frags = cellfrags.readlines()

numfrags = 664 # WARNING: Number of fragments is temporarily hardcoded! Would be provided by the code that generates p4frag.
numfatoms = 12 # WARNING: Number of atoms per fragment is temporarily hardcoded. Would be provided by the code that generates p4frag.
frg_separator = "--" # Fragment separator string.
fcounter = 1 # Counter of processed fragments.
lcounter = 0 # Counter of lines.
HeaderLine = True # Flag to indicate the first line of a fragment .xyz file.

for line in frags:
    lcounter += 1
    if line.startswith(frg_separator):
        fcounter += 1
        HeaderLine = True
    else:
        ffnidx = "f" + str(fcounter).zfill(len(str(numfrags))) + ".xyz"
        with open(ffnidx, "a") as frgxyz:
            if HeaderLine == True:
                frgxyz.write(str(numfatoms) + "\n"+ "Fragment " + str(fcounter) + "\n")
                HeaderLine = False
            if HeaderLine == False:
                frgxyz.write(line)
