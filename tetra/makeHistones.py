#!/usr/bin/env python
# Author:  Steven C. Howell
# Purpose: Prepare PDB for modeling
# Created: 04/24/2014
# $Id: makeHistones.py,v 1.1 2014-05-23 13:57:13 schowell Exp $
'''
This script creates a separate pdb for each chain of 'aa_pdb'
It also creates a sequence file for the residue sequence of that chain

'''
import sys
import os
import os.path as op
import subprocess
import logging
import sassie.sasmol.sasmol as sasmol
import numpy as np

class residue:
    def __init__(self, resid, resname):
        self.resid = resid
        self.resname = resname
    def __repr__(self):
        return repr((self.resid, self.resname))
        
aa_pdb = 'H2A_1kx5.pdb'
# aa_pdb = '../1zbb_tetra_uncombined.pdb'
# aa_pdb = '1zbb_original.pdb'

aa = sasmol.SasMol(0)
aa.read_pdb(aa_pdb)
chain_mols = []
errors = []


amino_acids = {'ALA': 'A',
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

dna = {'DG': 'G',
       'DA': 'A',
       'DT': 'T',
       'DC': 'C', 
       'GUA': 'G', 
       'ADE': 'A', 
       'THY': 'T', 
       'CYT': 'C'}

chains = ['C']


for chain in chains:
    if chain.lower() == chain:
        chain_name = 'chain_' + chain + '0'
    else:
        chain_name = 'chain_' + chain + '1'
    chain_name = '1kx5_' + chain_name
    basis_filter = "(chain[i] == '" + chain + "')"
    error, mask = aa.get_subset_mask(basis_filter)
    if error: print error
    chain_mol = sasmol.SasMol(0)
    error = aa.copy_molecule_using_mask(chain_mol, mask, 0)
    if error: print error
    chain_mol.write_pdb(chain_name+'.pdb', 0, 'w')
    chain_mols.append(chain_mol)

    # resids.sort()
    resA = 0
    res_min = np.min(chain_mol.resids())
    res_max = np.max(chain_mol.resids())
    print 'min resid:', res_min
    
    resA = 0
    residue_list = []
    
    ## create a sorted list of the residues
    for (i, resB) in enumerate(chain_mol.resid()):
        if resB != resA:
            # print 'chain_mol.resname()[i]:', chain_mol.resname()[i]
            # print 'chain_mol.resid()[i]:', chain_mol.resid()[i]
            residue_list.append(residue(resB, chain_mol.resname()[i]))
            resA = resB
        
    residue_sequence = sorted(residue_list, key=lambda residue: residue.resid)
    #with open(chain_name+'.txt', 'w') as outFile:
        #for res in residue_sequence:
            #outFile.write(str(res.resid) + '\t' + res.resname + '\n')

    with open(chain_name+'.seq', 'w') as outFile:          
        if chain_mol.moltypes() == ['protein']:
            for (i, res) in enumerate(residue_sequence):
                outFile.write(amino_acids[res.resname])
                if 0 == (i + 1) % 50:
                    outFile.write('\n')
        elif chain_mol.moltypes() == ['dna']:
            for (i, res) in enumerate(residue_sequence):
                outFile.write(dna[res.resname])
                if 0 == (i + 1) % 50:
                    outFile.write('\n')
        else:
            print 'ERROR, unexpected molecule type'
        
    #s for resB in resids:
    #s     if resB != resA + 1:
    #s         print 'missing residue/s, skipped chain', chain, 'btwn residues:', resA, resB
    #s     resA = resB
    print 'max resid:', res_max
    
    print 'finished chain', chain_name


print 'COMPLETE'
    
    
