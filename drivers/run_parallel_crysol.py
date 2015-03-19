#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: run parallel_crysol
# Created: 19 March 2015
#
# $Id: $
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sassie_1_na.drivers.parallel_crysol as par_crysol

class inputs():
    def __init__(self, parent = None):
        pass



if __name__ == '__main__':
    in_vars = inputs()
    in_vars.ncpu = 4
    in_vars.runname = 'test'
    in_vars.pdb = '3x167.pdb'
    in_vars.dcd = '%s/dna_mc/%s.dcd' % (in_vars.runname, in_vars.pdb[:-4])
    in_vars.dcd = '%s/dna_mc/%s.dcd' % (in_vars.runname, 'test')
    in_vars.driver = '/home/schowell/data/code/pylib/sassie_1_na/drivers/my_crysol_driver.py'
    in_vars.maxh = 21 # Tetramer_maxh ~= (0.2 x 325)/pi     
    in_vars.fib  = 18 # more does not slow crysol down noticeably
    in_vars.maxs = 0.2
    in_vars.numpoints = 26
    in_vars.sleep = 10
    in_vars.debug = True
    # in_vars.maxh = 3 # Tetramer_maxh ~= (0.2 x 325)/pi     
    # in_vars.fib  = 5 # more does not slow crysol down noticeably

    par_crysol.main(in_vars)
    
    print 'COMPLETE\n\m/ >.< \m/'