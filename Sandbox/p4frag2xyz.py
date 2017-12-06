#!/usr/bin/env python

p4frag = 'bzfrag.p4'

with open(p4frag) as cellfrags:
    frags = cellfrags.readlines()
    #print(frags)

frg_separator = '--'
i = 0

for line in frags:
    if line.startswith(frg_separator):
        i += 1
    else:
        fnindex = str(i).zfill(5)
        print(fnindex)
