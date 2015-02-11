#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Align several structures using an align basis
# Created: 6 February 2015
#
# $Id$
# 
'''
This script loads a pdb structure file of DNA, and creates a '*.patches' file
with the psfgen patches needed to use psfgen to create the structure.
After running this script, the patches can be pasted into a psfgen file.
'''

import sassie.sasmol.sasmol as sasmol
import sassie.tools.align2 as a2
import logging, sys
import numpy as np

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-g", "--goal",    help="goal pdb")
    parser.add_argument("-r", "--ref",     help="reference pdb")
    parser.add_argument("-m", "--move",    help="pdb/dcd to align")
    parser.add_argument("-o", "--out",     help="output dcd file")
    parser.add_argument("-s", "--segname", help="segname to match")
    parser.add_argument("-p", "--path",    help="output path")
    parser.add_argument("-mn", "--min",    help="minimun residue to match")
    parser.add_argument("-mx", "--max",    help="minimun residue to match")

    return parser.parse_args()

def main():
    aa_goal_pdb    = ARGS.goal
    aa_move_pdb    = ARGS.ref
    aa_move_file   = ARGS.move
    save_file      = ARGS.out
    match_segname  = ARGS.segname
    path           = ARGS.path
    match_segname  = ARGS.segname
    match_res_max  = ARGS.max
    match_res_min  = ARGS.min

    # create the SasMol objects
    aa_move  = sasmol.SasMol(0)
    sub_goal = sasmol.SasMol(0)
    sub_move = sasmol.SasMol(0)
    aa_goal  = sasmol.SasMol(0)
    aa_move  = sasmol.SasMol(0)

    aa_goal.read_pdb(aa_goal_pdb)
    aa_move.read_pdb(aa_move_pdb)

    if aa_move_file[-3:] == 'pdb':
        aa_move.read_pdb(aa_move_file)
        n_frames = aa_move.number_of_frames()
        in_type = 'pdb'
    elif aa_move_file[-3:] == 'dcd':
        dcd_file = aa_move.open_dcd_read(aa_move_file)
        n_frames = dcd_file[2]
        in_type = 'dcd'
    else:
        message = "\n~~~ ERROR, unknown input type ~~~\n"
        print_failure(message, txtOutput) 
        return

    dcd_out_file = None
    out_type = save_file[-3:]
    if 'dcd' == out_type:
        dcd_out_file = aa_move.open_dcd_write(path+save_file)

    basis_filter = "((segname[i] == '%s') and (name[i] == 'CA') and (resid[i] >= %s) and (resid[i] <= %s))" % (match_segname, match_res_min, match_res_max)

    error, goal_seg_mask = aa_goal.get_subset_mask(basis_filter)
    error, move_seg_mask = aa_move.get_subset_mask(basis_filter)

    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)

    com_sub_goal = sub_goal.calccom(0)                            # calculate the center of mass of the subset of m1
    sub_goal.center(0)                                           # center the m1 coordinates
    coor_sub_goal = sub_goal.coor()[0]                            # get the m1 centered coordinates

    for i in xrange(n_frames):
        if in_type == 'dcd':
            aa_move.read_dcd_step(dcd_file, i)
            aa_move.center(0)                                            # move m2 to be centered at the origin
            error, sub_move.coor = aa_move.get_coor_using_mask(0, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            com_sub_move = sub_move.calccom(0)                            # calculate the center of mass of the subset of m2
            sub_move.center(0)                                           # move the subset of m2 to be centered at the origin
            coor_sub_move = sub_move.coor[0]                              # get the new coordinates of the subset of m2
            aa_move.align(i,coor_sub_move,com_sub_move,coor_sub_goal,com_sub_goal) # align m2 using the transformation from sub_m2 to sub_m1
        elif in_type == 'pdb':
            aa_move.center(i)                                            # move m2 to be centered at the origin
            error, sub_move.coor = aa_move.get_coor_using_mask(i, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            com_sub_move = sub_move.calccom(0)                            # calculate the center of mass of the subset of m2
            sub_move.center(0)                                           # move the subset of m2 to be centered at the origin
            coor_sub_move = sub_move.coor[0]                              # get the new coordinates of the subset of m2
            aa_move.align(i,coor_sub_move,com_sub_move,coor_sub_goal,com_sub_goal) # align m2 using the transformation from sub_m2 to sub_m1

        a2.write_frame_to_file(aa_move,n_frames, i+1, path, save_file, out_type, dcd_out_file, in_type)

    if(intype == 'dcd'):
        aa_move.close_dcd_read(dcd_file[0])
    if(outtype == 'dcd'):
        aa_move.close_dcd_write(dcd_out_file)
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