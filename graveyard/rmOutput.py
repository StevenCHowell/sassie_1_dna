#!/usr/bin/python
#
# Author:  Steven C. Howell
# Purpose: remove debug files used frequently but not important
# Created: January 2014
#
# $Id$
#
import os

files = ('cg_dna.xyz', 'cg_dna_moves.dcd', 'cg_test.pdb', 'cg_test.dcd')
for fName in files:
    try:
        os.remove(fName)
        print 'removed: ' + fName
    except:
        print 'no such file: ' + fName
