#!/usr/bin/python
# $Id: rmOutput.py,v 1.1 2013-11-20 16:29:27 schowell Exp $
import os

files = ('cg_dna.xyz','cg_dna_moves.dcd','cg_test.pdb','cg_test.dcd')
for fName in files:
    try:
        os.remove(fName)
        print 'removed: ' + fName
    except :
        print 'no such file: ' + fName
