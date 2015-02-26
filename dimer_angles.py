#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Plot the twist angle for the dinucleosome
# Created: 26 February 2015
#
# $Id$
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sys
# import os
import os.path as op

import logging
LOGGER = logging.getLogger(__name__) #add module name manually

import sassie_1_na.fit_cylinder as fit_cyl


class MainError(Exception):
    pass


def main():
    NotImplemented


if __name__ == '__main__':

    if '-v' in sys.argv:
        logging.basicConfig(filename='%s.log' %__file__[:-3], 
                            level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    main()    
    
    
    
    
    print '\m/ >.< \m/'