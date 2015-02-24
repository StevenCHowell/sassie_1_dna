#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Fit a cylider to a set of coordinates
# Created: 23 February 2015
#
# $Id$
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sys
import os
import os.path as op
import argparse #http://docs.python.org/dev/library/argparse.html
import numpy as np

import logging
LOGGER = logging.getLogger(__name__) #add module name manually


class MainError(Exception):
    pass

def get_c1p_coors():
    NotImplemented
    
def calc_distances(coor, R, X0, Y0, Vx, Vy):
    Z0 = 0
    Vz = 1
    length = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    Ux = Vx/length
    Uy = Vy/length
    Uz = Vz/length
    A = Ux*(coor-)

def main():
    NotImplemented


def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        prog='',
        usage='',
        description = 'description here',
        epilog = 'no epilog found'
        )

    parser.add_argument("positional", help='')
    parser.add_argument("--option", metavar='', help='')
    return parser.parse_args()


if __name__ == '__main__':

    if '-v' in sys.argv:
        logging.basicConfig(level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

   # make ARGS global
    ARGS = parse()
    
    
    main()