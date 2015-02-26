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

    dcd_file = 'dimer_ta.dcd'
    dcd_dir  = '/home/schowell/data/myData/sassieRuns/dimer/flex25/run2'
    dcd_path = op.join(dcd_dir, dcd_file)

    pdb_file = 'dimer_half_trun.pdb'
    pdb_dir  = '/home/schowell/data/myData/sassieRuns/dimer/flex25/run2'
    pdb_path = op.join(pdb_dir, pdb_path)
    
    # load the pdb
    # create the ncp filters
    
    for frame in xrange(20):
        # load the coordinates from the dcd frame
        
        ncp1_origin, ncp1_axes = fit_cyl.get_ncp_origin_and_axes(ncp_c1p_filter, dyad_dna_resids, dyad_dna_id, ncp)
        
        ncp2_origin, ncp2_axes = fit_cyl.get_ncp_origin_and_axes(ncp_c1p_filter, dyad_dna_resids, dyad_dna_id, ncp)

        #get angles between NCPs

    #get the X2 from the file

    #plot the X2 vs angles

    
    main()    
    
    
    print '\m/ >.< \m/'