#!/usr/bin/env python
# Author:  Steven C. Howell
# Purpose: Replace DNA sequence with another sequence
# Created: 04/24/2014
# $Id: pdb2psfgen.py,v 1.2 2014-05-23 14:15:03 schowell Exp $
'''
This script loads a pdb structure file of DNA, and creates a '*.patches' file
with the psfgen patches needed to use psfgen to create the structure.
After running this script, the patches can be pasted into a psfgen file.
'''

import sassie.sasmol.sasmol as sasmol
import sys, string
import time

allAtomPDB = sys.argv[1]
m1 = sasmol.SasMol(0)
m1.read_pdb(allAtomPDB)

chain1 = 'J'
chain2 = 'L'
    
names = m1.resname()
ids = m1.resid()
c = m1.chain()

pdbStrings = ['.pdb','.PDB']
psfgenFile = sys.argv[1]
#psfgenFile = sys.argv[2]
for ext, item in enumerate(pdbStrings):
    psfgenFile = psfgenFile.replace(item, '.patches')

outfile = open(psfgenFile, 'w')      # open the file
timestr = time.strftime("created on %d %B %Y by 'pdb2psfgen.py'")
outfile.write(timestr)

pyr = ['DC', 'DT', 'CYT', 'THY']
pur = ['DA', 'DG', 'ADE', 'GUA']
pyrStr = 'patch DEO1 '
purStr = 'patch DEO2 '

n = 0
for (j, i) in enumerate(ids):
    # only want this to happend once for each residue
    if n != i:
        n = i
        # print 'adding line %d' % i
        
        if c[j] in chain1:
            dna = 'dna1:%d\n' %i
        elif c[j] in chain2:
            dna = 'dna2:%d\n' %i
        else:
            print 'ERROR!!! residue from unknown dna chain (protein ?)'
            dna = 'protein:%d\n' % i
            
        if names[j] in pyr:
            outfile.write(pyrStr + dna)
            # print pyrStr + dna
        elif names[j] in pur:
            outfile.write(purStr + dna)
            # print purStr + dna
        else:
            print '\n\n ERROR!!! skipped unknown resname'

outfile.close()            

print 'COMPLETE'