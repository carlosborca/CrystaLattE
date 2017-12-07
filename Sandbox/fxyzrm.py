#!/usr/bin/env python

import os
import re

numfatoms = 12 # WARNING: Number of atoms per fragment hardcoded!
               # WARNING: What if there are two types of molecules?

directory = os.getcwd()
files = os.listdir(directory) 

print("")
print("Detecting fragments with incomplete molecules:")
print("")

incmolcounter = 0

for file in files:
    if re.match('^f[0-9]+.xyz$', file):
        with open(file, 'r') as f:
            lines = f.readlines()
            if len(lines) != numfatoms + 2:
                print("Expected " + str(numfatoms) + " atoms. Found only " + str(len(lines) - 2) + ". Removing: " + file)
                os.remove(os.path.join(directory,file))
                incmolcounter += 1

print ("Removed %s fragments." % incmolcounter)
