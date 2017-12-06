#!/usr/bin/env python

p4frag = "bzfrag.p4"

with open(p4frag) as cellfrags:
    frags = cellfrags.readlines()
    #print(frags)

numfrags = 664 # WARNING: Number of fragments is hardcoded temporarily! Would be provided by the code that generates p4frag.
numfatoms = 12 # WARNING: Number of atoms per fragment is hardcoded temporarily. Would be provided by the code that generates p4frag.
frg_separator = "--"
fcounter = 1
lcounter = 0
HeaderLine = True

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
