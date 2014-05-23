#!/usr/bin/python
import sassie.sasmol.sasmol as sasmol
import sys, string
import time

allAtomPDB = sys.argv[1]
m1 = sasmol.SasMol(0)
m1.read_pdb(allAtomPDB)


    
names = m1.resname()
ids = m1.resid()
c = m1.chain()

pdbStrings = ['.pdb','.PDB']
psfgenFile = sys.argv[1]
#psfgenFile = sys.argv[2]
for ext, item in enumerate(pdbStrings):
    psfgenFile = psfgenFile.replace(item, '.patchs')

outfile = open(psfgenFile, 'w')      # open the file
timestr = time.strftime("created on %d %B %Y by 'pdb2psfgen.py'")
outfile.write(timestr)

pyr = ['DC', 'DT', 'CYT', 'THY']
pur = ['DA', 'DG', 'ADE', 'GUA']
pyrStr = 'patch DEO1 '
purStr = 'patch DEO2 '

n = 0
j = 0
for j in xrange(ids.size):
    i = ids[j]
    if n != i:
        n = i
        # print 'adding line %d' % i
        
        if c[j] in ['A']:
            dna = 'dna1:%d\n' %i
        elif c[j] in ['B']:
            dna = 'dna2:%d\n' %i
        else:
            print 'ERROR!!! residue from unknown dna chain'
            dna = 'dna?:%d' % i
            
        if names[j] in pyr:
            outfile.write(pyrStr + dna)
            # print pyrStr + dna
        elif names[j] in pur:
            outfile.write(purStr + dna)
            # print purStr + dna
        else:
            print '\n\n ERROR!!! skipped unknown resname'

outfile.close()            
