#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Prepare PDB for modeling
# Created: 24 April 2014
#
# $Id$
#
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


def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        # prog='',
        # usage='',
        description='test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-p", "--pdb", help="all atom pdb file")
    parser.add_argument(
        "-s", "--segnames", nargs='+', help="segnames to extract")

    return parser.parse_args()


def main():
    # aa_pdb = '../1zbb_tetra_uncombined.pdb'
    # aa_pdb = '1zbb_original.pdb'

    aa = sasmol.SasMol(0)
    aa.read_pdb(ARGS.pdb)
    segname_mols = []
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
                   'VAL': 'V',
                   'HSE': 'H'}

    dna = {'G': 'G',
           'A': 'A',
           'T': 'T',
           'C': 'C',
           'DG': 'G',
           'DA': 'A',
           'DT': 'T',
           'DC': 'C',
           'GUA': 'G',
           'ADE': 'A',
           'THY': 'T',
           'CYT': 'C'}

    # segnames = ['I','J']
    # print 'segnames =', segnames
    print 'ARGS.segnames =', ARGS.segnames

    for segname in ARGS.segnames:
        if segname.lower() == segname:
            segname_name = '_seg_' + segname + '0'
        else:
            segname_name = '_seg_' + segname + '1'
        segname_name = ARGS.pdb[:-4] + segname_name
        basis_filter = "(segname[i] == '" + segname + "')"
        error, mask = aa.get_subset_mask(basis_filter)
        if error:
            print error
        segname_mol = sasmol.SasMol(0)
        error = aa.copy_molecule_using_mask(segname_mol, mask, 0)
        if error:
            print error
        segname_mol.write_pdb(segname_name + '.pdb', 0, 'w')
        segname_mols.append(segname_mol)

        # resids.sort()
        resA = 0
        res_min = np.min(segname_mol.resids())
        res_max = np.max(segname_mol.resids())
        print 'min resid:', res_min

        resA = 0
        residue_list = []

        # create a sorted list of the residues
        for (i, resB) in enumerate(segname_mol.resid()):
            if resB != resA:
                # print 'segname_mol.resname()[i]:', segname_mol.resname()[i]
                # print 'segname_mol.resid()[i]:', segname_mol.resid()[i]
                residue_list.append(residue(resB, segname_mol.resname()[i]))
                resA = resB

        residue_sequence = sorted(
            residue_list, key=lambda residue: residue.resid)
        # with open(segname_name+'.txt', 'w') as outFile:
        # for res in residue_sequence:
        #outFile.write(str(res.resid) + '\t' + res.resname + '\n')
        if 'rna' in segname_mol.moltypes():
            segname_mol.moltypes().remove('rna')
            print "removed 'rna' from moltypes"
            if len(segname_mol.moltypes()) == 0:
                segname_mol.moltypes().append('dna')
                print "appended 'dna' to moltypes"

        with open(segname_name + '.seq', 'w') as outFile:
            print outFile.closed
            if segname_mol.moltypes() == ['protein']:
                for (i, res) in enumerate(residue_sequence):
                    outFile.write(amino_acids[res.resname])
                    if 0 == (i + 1) % 50:
                        outFile.write('\n')
            elif segname_mol.moltypes() == ['dna']:
                for (i, res) in enumerate(residue_sequence):
                    outFile.write(dna[res.resname])
                    # print 'printed', dna[res.resname], 'to', segname_name,
                    # '.seq'
                    if 0 == (i + 1) % 50:
                        outFile.write('\n')
            else:
                print 'ERROR, unexpected molecule type'

        # s for resB in resids:
        # s     if resB != resA + 1:
        # s         print 'missing residue/s, skipped segname', segname, 'btwn residues:', resA, resB
        # s     resA = resB
        print 'max resid:', res_max

        print 'finished segname', segname_name
        print outFile.closed

    print 'COMPLETE'

if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()
    main()

    print '!!! \m/ >.< \m/ !!!'
