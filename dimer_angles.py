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
from numpy.core.umath_tests import inner1d
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

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
    
    aligned = True
    
    all_ncp1_origins = []
    all_ncp1_axes = []
    all_ncp2_origins = []
    all_ncp2_axes = []
    
    import time
    tic = time.time()
    n_frames = 1
    
    ncp2_opt_params = None
    
    for frame in xrange(n_frames):
        # load the coordinates from the dcd frame
        dimer.read_dcd_step(dimer_dcdfile, frame)
        
        if not aligned or frame is 0:
            ncp1_origin, ncp1_axes, ncp1_opt_params = fit_cylinder.get_ncp_origin_and_axes(
                ncp1_c1p_mask, ncp1_dyad_resids, dna_ids, dimer, dna_id_type='segname', debug=False)
            all_ncp1_origins.append(ncp1_origin)
            all_ncp1_axes.append(ncp1_axes)
            toc = time.time() - tic
            print 'fitting ncp1 took %0.3f s' %toc    
        
        try:
            ncp2_origin, ncp2_axes, ncp2_opt_params = fit_cylinder.get_ncp_origin_and_axes(
                ncp2_c1p_mask, ncp2_dyad_resids, dna_ids, dimer, ncp2_opt_params, 'segname', debug=False)
        except RuntimeError:
            print 'ERROR: curve_fit failed'
            ncp2_origin = all_ncp2_origins[-1] *  0
            ncp2_axes = all_ncp2_axes[-1] * 0
        all_ncp2_origins.append(ncp2_origin)
        all_ncp2_axes.append(ncp2_axes)
    toc = time.time() - tic
    print 'fitting both NCPs for %d iterations took %0.3f s' % (n_frames, toc)

    all_ncp2_axes = np.array(all_ncp2_axes)
    all_ncp2_origins = np.array(all_ncp2_origins)

    # get the minimum distance from the origin and the axis
    ncp1_axis_vectors = vector_from_line(dimer.coor(), ncp1_origin, ncp1_axes[2])
    ncp1_axis_magnitude = np.sqrt(inner1d(ncp1_axis_vectors, ncp1_axis_vectors))
    ncp1_origin_vectors = dimer.coor() - ncp1_origin
    ncp1_origin_magnitude = np.sqrt(inner1d(ncp1_origin_vectors, ncp1_origin_vectors))
    
    ncp2_axis_vectors = vector_from_line(dimer.coor(), ncp2_origin, ncp2_axes[2])
    ncp2_axis_magnitude = np.sqrt(inner1d(ncp2_axis_vectors, ncp2_axis_vectors))
    ncp2_origin_vectors = dimer.coor() - ncp2_origin
    ncp2_origin_magnitude = np.sqrt(inner1d(ncp2_origin_vectors, ncp2_origin_vectors))

    all_ncp1_axes = np.copy(all_ncp2_axes)
    all_ncp1_origins = np.copy(all_ncp2_origins)
    all_ncp1_axes[:] = ncp1_axes
    all_ncp1_origins[:] = all_ncp1_origins

    ## get angles between NCPs
    # bending
    phi = np.arccos(inner1d(all_ncp1_axes[:,0], all_ncp2_axes[:,0]))
    # twist (Victor's group from NIH uses rho = psi/2)
    psi = np.arccos(inner1d(all_ncp1_axes[:,2], all_ncp2_axes[:,2]))

    plt.scatter(phi,psi)
    plt.xlabel('phi: bend')
    plt.ylabel('psi: twist')
    plt.show()
    
    #get the X2 from the file

    #plot the X2 vs angles

    
    main()    
    
    
    print '\m/ >.< \m/'