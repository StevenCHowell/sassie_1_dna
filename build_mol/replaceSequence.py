#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Replace DNA sequence with another sequence
# Created: 24 April 2014
#
# $Id$
#
'''
This script loads a pdb structure file of DNA, replaces the DNA sequence,
then saves the result as a new pdb.  The pdb input should just be the DNA
backbone (this can be achieved by loading a pdb into VMD, then saving just
the backbone of the DNA into a new PDB). It is important to check that the
same number of residues are in the two inputs.
'''



import sys
import os
import os.path as op
import subprocess
import logging
import sassie.sasmol.sasmol as sasmol
import numpy as np

chain1_out = 'chainPDB/c11_chain1_backbone.pdb'
chain1_pdb = 'chainPDB/tetramer_chain1_backbone.pdb'
chain1 = sasmol.SasMol(0)
chain1.read_pdb(chain1_pdb)

chain2_out = 'chainPDB/c11_chain2_backbone.pdb'
chain2_pdb = 'chainPDB/tetramer_chain2_backbone.pdb'
chain2 = sasmol.SasMol(0)
chain2.read_pdb(chain2_pdb)

sequenceFile = 'sequences/c11.seq'
with open(sequenceFile) as f:
    lines = f.read().splitlines()
sequence = lines[0]
reverse = sequence[::-1]

amino_acids_pdb2seq = {'ALA': 'A',
                       'ARG': 'R',
                       'ASN': 'N',
                       'ASP': 'D',
                       'CYS': 'C',
                       'GLU': 'E',
                       'GLN': 'Q',
                       'GLY': 'G',
                       'HIS': 'H',
                       'ILE': 'I',
                       'LEU': 'L',
                       'LYS': 'K',
                       'MET': 'M',
                       'PHE': 'F',
                       'PRO': 'P',
                       'SER': 'S',
                       'THR': 'T',
                       'TRP': 'W',
                       'TYR': 'Y',
                       'VAL': 'V'}

dna_pdb2seq = {'DG' : 'G',
               'DA' : 'A',
               'DT' : 'T',
               'DC' : 'C', 
               'GUA': 'G', 
               'ADE': 'A', 
               'THY': 'T', 
               'CYT': 'C'}

dna_seq2pdb = {'G': 'GUA', 
               'A': 'ADE', 
               'T': 'THY', 
               'C': 'CYT'}

dna_compliment = {'G': 'C', 
                  'A': 'T', 
                  'T': 'A', 
                  'C': 'G'}


for (i, atom) in enumerate(chain1.resid()):
    # print 'changing residue ', aa.resname()[i], ' from ', aa.res
    # print 'atom: ', atom
    print atom, '- chainging aa.resname()[', i, '] from: ', chain1.resname()[i], 'to: ', dna_seq2pdb[sequence[atom-1]]
    # print 'sequence[atom-1]: ', sequence[atom-1]
    # print 'dna_seq2pdb[sequence[atom-1]]: ', dna_seq2pdb[sequence[atom-1]]
    chain1.resname()[i] = dna_seq2pdb[sequence[atom-1]]
    
for (i, atom) in enumerate(chain2.resid()):
    print 'atom: ', atom
    compliment = dna_compliment[reverse[atom-1]]
    print 'compliment: ', compliment

    chain2.resname()[i] = dna_seq2pdb[compliment]
    
chain1.setResname(chain1.resname())
chain2.setResname(chain2.resname())

chain1.write_pdb(chain1_out, 0, 'w')
chain2.write_pdb(chain2_out, 0, 'w')


print 'COMPLETE'
    
    
