#!/usr/bin/env python

import Read_CIF

args = ['', '-i', 'Crystal_Benzene-138K.cif', '-o',  'Crystal_Benzene-138K.xyz', '-b', '5', '5', '5']

Read_CIF.main(args)
