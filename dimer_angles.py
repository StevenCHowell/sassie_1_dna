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

import sassie.sasmol.sasmol as sasmol
import sassie_1_na.util.fit_cylinder as fit_cylinder
import sassie_1_na.util.basis_to_python as basis_to_python
import numpy as np

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
    pdb_dir  = '/home/schowell/data/myData/sassieRuns/dimer/flex25'
    pdb_path = op.join(pdb_dir, pdb_file)
    
    # load the pdb and open the pdb for reading 
    dimer = sasmol.SasMol(0)
    dimer.read_pdb(pdb_path)
    dimer_dcdfile = dimer.open_dcd_read(dcd_path)
    
    # create the ncps
    ncp1_vmd_basis = ("((segname DNA1 and resid >= 15  and resid <= 161) or "
                      " (segname DNA2 and resid >= 193 and resid <= 339) ) and "
                      "name C1'")
    ncp1_basis = basis_to_python.parse_basis(ncp1_vmd_basis)
    error, ncp1_c1p_mask = dimer.get_subset_mask(ncp1_basis)

    ncp2_vmd_basis = ("((segname DNA1 and resid >= 182 and resid <= 328) or "
                      " (segname DNA2 and resid >= 26  and resid <= 172) )and "
                      "name C1'")
    ncp2_basis = basis_to_python.parse_basis(ncp2_vmd_basis)
    error, ncp2_c1p_mask = dimer.get_subset_mask(ncp2_basis)

    ncp1_dyad_resids = [88, 266]
    ncp2_dyad_resids = [255, 99]
    dna_ids = ['DNA1', 'DNA2']
    
    all_ncp1_origins = []
    all_ncp1_axes = []
    all_ncp2_origins = []
    all_ncp2_axes = []
    for frame in xrange(20):
        # load the coordinates from the dcd frame
        dimer.read_dcd_step(dimer_dcdfile, frame)
        ncp1_origin, ncp1_axes = fit_cylinder.get_ncp_origin_and_axes(ncp1_c1p_mask, ncp1_dyad_resids, dna_ids, dimer, 'segname')
        all_ncp1_origins.append(ncp1_origin)
        all_ncp1_axes.append(ncp1_axes)
        
        ncp2_origin, ncp2_axes = fit_cylinder.get_ncp_origin_and_axes(ncp2_c1p_mask, ncp2_dyad_resids, dna_ids, dimer, 'segname')
        all_ncp2_origins.append(ncp2_origin)
        all_ncp2_axes.append(ncp2_axes)

    all_ncp1_axes = np.array(all_ncp1_axes)
    all_ncp1_origins = np.array(all_ncp1_axes)
    all_ncp2_axes = np.array(all_ncp2_axes)
    all_ncp2_origins = np.array(all_ncp2_axes)
    #get angles between NCPs

    #get the X2 from the file

    #plot the X2 vs angles

    
    main()    
    
    
    print '\m/ >.< \m/'