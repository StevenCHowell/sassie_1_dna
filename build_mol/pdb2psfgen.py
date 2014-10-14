#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Replace DNA sequence with another sequence
# Created: 24 April 2014
#
# $Id$
# 
'''
This script loads a pdb structure file of DNA, and creates a '*.patches' file
with the psfgen patches needed to use psfgen to create the structure.
After running this script, the patches can be pasted into a psfgen file.
'''

import sassie.sasmol.sasmol as sasmol
import sys, string
import time
import logging

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-p", "--pdb", help="all atom pdb file")
    parser.add_argument("-c", "--chains", nargs='+',  
        help="chains to extract")

    return parser.parse_args()

def main():
    m1 = sasmol.SasMol(0)
    m1.read_pdb(ARGS.pdb)
    
    print ARGS.chains    
    
    chain1 = ARGS.chains[0]
    chain2 = ARGS.chains[1]
        
    print 'chain 1: ', chain1
    print 'chain 2: ', chain2
    
    names = m1.resname()
    ids = m1.resid()
    c = m1.chain()
    
    psfgenFile = ARGS.pdb[:-4] + '_patches.txt'
    
    outfile = open(psfgenFile, 'w')      # open the file
    timestr = time.strftime("# created on %d %B %Y by 'pdb2psfgen.py'\n")
    outfile.write(timestr)
    outfile.write('# dna1: chain ' + chain1 + '\n')
    outfile.write('# dna2: chain ' + chain2 + '\n')
    
    pyr = ['C','T','DC', 'DT', 'CYT', 'THY']
    pur = ['A','G','DA', 'DG', 'ADE', 'GUA']
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
                print 'Skipping residue from unspecified chain: ', c[j]
                break
                #s dna = 'protein:%d\n' % i
                
            if names[j] in pyr:
                outfile.write(pyrStr + dna)
                # print pyrStr + dna
            elif names[j] in pur:
                outfile.write(purStr + dna)
                # print purStr + dna
            else:
                print 'ERROR!!! unknown resname in specified chain: ', names[j]
                print '\n'
    
    outfile.close()            
    
    print 'COMPLETE \m/ >.< \m/'
    
if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()
    main()
