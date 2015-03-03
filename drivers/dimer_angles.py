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
import sassie_1_na.util.geometry as geometry
import sassie_1_na.util.basis_to_python as basis_to_python
import numpy as np
from numpy.core.umath_tests import inner1d
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

class MainError(Exception):
    pass


def main(dimer, dimer_dcdfile):
    # create the ncps
    ncp1_vmd_basis = ("((segname DNA1 and resid >= 15  and resid <= 161) or "
                      " (segname DNA2 and resid >= 193 and resid <= 339) ) and "
                      "name C1'")
    ncp1_basis = basis_to_python.parse_basis(ncp1_vmd_basis)
    error, ncp1_c1p_mask = dimer.get_subset_mask(ncp1_basis)

    ncp2_vmd_basis = ("((segname DNA1 and resid >= 182 and resid <= 328) or "
                      " (segname DNA2 and resid >= 26  and resid <= 172) )and "
                      "backbone'")
    ncp2_basis = basis_to_python.parse_basis(ncp2_vmd_basis)
    error, ncp2_c1p_mask = dimer.get_subset_mask(ncp2_basis)

    ncp1_dyad_resids = [88, 266]
    ncp2_dyad_resids = [255, 99]
    dna_ids = ['DNA1', 'DNA2']
    
    aligned = True
    
    import time
    n_frames = 1000
    
    tic = time.time()
    all_ncp1_origins, all_ncp2_origins, all_ncp1_axes, all_ncp2_axes = fit_method(
        n_frames, dimer_dcdfile, dimer, aligned, ncp1_c1p_mask, ncp1_dyad_resids, 
        dna_ids, ncp2_dyad_resids, ncp2_c1p_mask)
    toc = time.time() - tic
    print 'fitting both NCPs for %d iterations took %0.3f s' % (n_frames, toc)

    ## get angles between NCPs
    # bending
    phi_fit_r = np.arccos(inner1d(all_ncp1_axes[:,0], all_ncp2_axes[:,0]))
    phi_fit_d = phi_fit_r * 180 / np.pi
    # twist (Victor's group from NIH uses rho = psi/2)
    psi_fit_r = np.arccos(inner1d(all_ncp1_axes[:,2], all_ncp2_axes[:,2]))
    psi_fit_d = psi_fit_r * 180 / np.pi
    

    # tic = time.time()
    # all_ncp1_axes, all_ncp2_axes = mask_method(ncp1_c1p_mask, ncp1_dyad_resids, 
                                               # ncp2_c1p_mask, ncp2_dyad_resids, 
                                               # aligned, dimer_dcdfile, dna_ids, 
                                               # dimer, n_frames)
    # toc = time.time() - tic
    # print 'getting NCP axes using masks for %d iterations took %0.3f s' % (n_frames, toc)

    ## get angles between NCPs
    # bending
    phi_mask_r = np.arccos(inner1d(all_ncp1_axes[:,0], all_ncp2_axes[:,0]))
    phi_mask_d = phi_mask_r * 180 / np.pi
    # twist (Victor's group from NIH uses rho = psi/2)
    psi_mask_r = np.arccos(inner1d(all_ncp1_axes[:,2], all_ncp2_axes[:,2]))
    psi_mask_d = psi_mask_r * 180 / np.pi

    # plt.scatter(phi_fit_d,psi_fit_d, 'ro')
    # plt.scatter(phi_mask_d,psi_mask_d, 'b^')
    # plt.xlabel('phi: bend')
    # plt.ylabel('psi: twist')
    # plt.show()
    
    #get the X2 from the file
    
    #plot the X2 vs anglesdef main():    
    
    return

def mask_method(ncp1_c1p_mask, ncp1_dyad_resids, ncp2_c1p_mask, ncp2_dyad_resids, 
                aligned, dimer_dcdfile, dna_ids, dimer, n_frames):

    (i_min_ncp1_x, i_min_ncp1_y, i_min_ncp1_z, i_min_ncp1_origin, 
     i_min_ncp2_x, i_min_ncp2_y, i_min_ncp2_z, i_min_ncp2_origin) = get_closest_atoms(
         dimer_dcdfile, dimer, aligned, ncp1_c1p_mask, ncp1_dyad_resids, 
         dna_ids, ncp2_dyad_resids, ncp2_c1p_mask)    
    all_ncp1_origins = []
    all_ncp1_axes = []
    all_ncp2_origins = []
    all_ncp2_axes = []
    
    for frame in xrange(n_frames):
        # load the coordinates from the dcd frame
        dimer.read_dcd_step(dimer_dcdfile, frame)
        
        if not aligned or frame is 0:
            ncp1_x_coor = dimer.coor()[0,i_min_ncp1_x]
            ncp1_y_coor = dimer.coor()[0,i_min_ncp1_y]
            # ncp1_z_coor = dimer.coor()[0,i_min_ncp1_z]
            ncp1_origin = dimer.coor()[0,i_min_ncp1_origin]
            ncp1_x, ncp1_y, ncp1_z = geometry.get_axes_from_points(ncp1_origin, ncp1_x_coor, ncp1_y_coor)
            ncp1_axes = np.array([ncp1_x, ncp1_y, ncp1_z])  
        
        ncp2_x_coor = dimer.coor()[0,i_min_ncp2_x]
        ncp2_y_coor = dimer.coor()[0,i_min_ncp2_y]
        # ncp2_z_coor = dimer.coor()[0,i_min_ncp2_z]
        ncp2_origin = dimer.coor()[0,i_min_ncp2_origin]
        ncp2_x, ncp2_y, ncp2_z = geometry.get_axes_from_points(ncp2_origin, ncp2_x_coor, ncp2_y_coor)
        all_ncp2_origins.append(ncp2_origin)
        all_ncp2_axes.append(np.array([ncp2_x, ncp2_y, ncp2_z]))        

    all_ncp2_axes = np.array(all_ncp2_axes)
    all_ncp2_origins = np.array(all_ncp2_origins)

    if aligned:
        all_ncp1_axes = np.copy(all_ncp2_axes)
        all_ncp1_origins = np.copy(all_ncp2_origins)
        all_ncp1_axes[:] = ncp1_axes
        all_ncp1_origins[:] = ncp1_origin
    else:
        all_ncp1_axes = np.array(all_ncp1_axes)
        all_ncp1_origins = np.array(all_ncp1_origins)

    return all_ncp1_axes, all_ncp2_axes
    
def get_closest_atoms(dimer_dcdfile, dimer, aligned, ncp1_c1p_mask, ncp1_dyad_resids, 
        dna_ids, ncp2_dyad_resids, ncp2_c1p_mask):
    n_frames = 1
    
    all_ncp1_origins, all_ncp2_origins, all_ncp1_axes, all_ncp2_axes = fit_method(
        n_frames, dimer_dcdfile, dimer, aligned, ncp1_c1p_mask, ncp1_dyad_resids, 
        dna_ids, ncp2_dyad_resids, ncp2_c1p_mask, True)
    
    ncp1_origin = all_ncp1_origins[0]
    ncp1_axes = all_ncp1_axes[0]
    ncp2_origin = all_ncp2_origins[0]
    ncp2_axes = all_ncp2_axes[0]
    # get the minimum distance from the x and z axis (y obtiained by r-h-r)
    ncp1_x_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp1_origin, ncp1_axes[0])
    ncp1_x_axis_magnitudes = np.sqrt(inner1d(ncp1_x_axis_vectors, ncp1_x_axis_vectors))
    print_min(ncp1_x_axis_magnitudes, dimer, 'ncp1_x_axis')
    i_min_ncp1_x = ncp1_x_axis_magnitudes.argmin()
    ncp1_y_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp1_origin, ncp1_axes[1])
    ncp1_y_axis_magnitudes = np.sqrt(inner1d(ncp1_y_axis_vectors, ncp1_y_axis_vectors))
    print_min(ncp1_y_axis_magnitudes, dimer, 'ncp1_y_axis')
    i_min_ncp1_y = ncp1_y_axis_magnitudes.argmin()
    ncp1_z_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp1_origin, ncp1_axes[2])
    ncp1_z_axis_magnitudes = np.sqrt(inner1d(ncp1_z_axis_vectors, ncp1_z_axis_vectors))
    print_min(ncp1_z_axis_magnitudes, dimer, 'ncp1_z_axis')
    i_min_ncp1_z = ncp1_z_axis_magnitudes.argmin()
    
    ncp2_x_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp2_origin, ncp2_axes[0])
    ncp2_x_axis_magnitudes = np.sqrt(inner1d(ncp2_x_axis_vectors, ncp2_x_axis_vectors))
    print_min(ncp2_x_axis_magnitudes, dimer, 'ncp2_x_axis')
    i_min_ncp2_x = ncp2_x_axis_magnitudes.argmin()
    ncp2_y_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp2_origin, ncp2_axes[1])
    ncp2_y_axis_magnitudes = np.sqrt(inner1d(ncp2_y_axis_vectors, ncp2_y_axis_vectors))
    print_min(ncp2_y_axis_magnitudes, dimer, 'ncp2_y_axis')
    i_min_ncp2_y = ncp2_y_axis_magnitudes.argmin()
    ncp2_z_axis_vectors = geometry.vector_from_line(dimer.coor()[0], ncp2_origin, ncp2_axes[2])
    ncp2_z_axis_magnitudes = np.sqrt(inner1d(ncp2_z_axis_vectors, ncp2_z_axis_vectors))
    print_min(ncp2_z_axis_magnitudes, dimer, 'ncp2_z_axis')
    i_min_ncp2_z = ncp2_z_axis_magnitudes.argmin()
    
    # get the mimimum distance from the origin
    ncp1_origin_vectors = dimer.coor() - ncp1_origin
    ncp1_origin_magnitudes = np.sqrt(inner1d(ncp1_origin_vectors, ncp1_origin_vectors))
    print_min(ncp1_origin_magnitudes, dimer, 'ncp1_origin')
    i_min_ncp1_origin = ncp1_origin_magnitudes.argmin()
    
    ncp2_origin_vectors = dimer.coor()[0] - ncp2_origin
    ncp2_origin_magnitudes = np.sqrt(inner1d(ncp2_origin_vectors, ncp2_origin_vectors))
    print_min(ncp2_origin_magnitudes, dimer, 'ncp2_origin')
    i_min_ncp2_origin = ncp2_origin_magnitudes.argmin()
    return (i_min_ncp1_x, i_min_ncp1_y, i_min_ncp1_z, i_min_ncp1_origin, 
            i_min_ncp2_x, i_min_ncp2_y, i_min_ncp2_z, i_min_ncp2_origin)

def fit_method(n_frames, dimer_dcdfile, dimer, aligned, ncp1_c1p_mask, 
               ncp1_dyad_resids, dna_ids, ncp2_dyad_resids, ncp2_c1p_mask, debug=False):
    all_ncp1_origins = []
    all_ncp1_axes = []
    all_ncp2_origins = []
    all_ncp2_axes = []
    
    ncp1_dyad_filter = '( segname[i] == "%s" and resid[i] == %d ) or ( segname[i] == "%s" and resid[i] == %d )' % (dna_ids[0], ncp1_dyad_resids[0], dna_ids[1], ncp1_dyad_resids[1])
    ncp2_dyad_filter = '( segname[i] == "%s" and resid[i] == %d ) or ( segname[i] == "%s" and resid[i] == %d )' % (dna_ids[0], ncp2_dyad_resids[0], dna_ids[1], ncp2_dyad_resids[1])
    error, ncp1_dyad_mask = dimer.get_subset_mask(ncp1_dyad_filter)
    error, ncp2_dyad_mask = dimer.get_subset_mask(ncp2_dyad_filter)
    
    ncp2_opt_params = ncp1_dyad_mol = ncp2_dyad_mol = None

    for frame in xrange(n_frames):
        # load the coordinates from the dcd frame
        dimer.read_dcd_step(dimer_dcdfile, frame)
        coor = dimer.coor()
        coor[0] = geometry.transform_coor(coor[0], np.array([1,0,0]), np.array([0,0,0]))
        dimer.setCoor(coor)
        
        if not aligned or frame is 0:
            ncp1_origin, ncp1_axes, ncp1_opt_params, ncp1_dyad_mol = geometry.get_ncp_origin_and_axes(
                ncp1_c1p_mask, ncp1_dyad_mask, dna_ids, dimer, None, 'segname', ncp1_dyad_mol, debug)
            if not aligned:
                all_ncp1_origins.append(ncp1_origin)
                all_ncp1_axes.append(ncp1_axes)
        
        try:
            ncp2_origin, ncp2_axes, ncp2_opt_params, ncp2_dyad_mol = geometry.get_ncp_origin_and_axes(
                ncp2_c1p_mask, ncp2_dyad_mask, dna_ids, dimer, ncp2_opt_params, 'segname', ncp2_dyad_mol, debug)
        except RuntimeError:
            try:
                ncp2_origin, ncp2_axes, ncp2_opt_params, ncp2_dyad_mol = geometry.get_ncp_origin_and_axes(
                    ncp2_c1p_mask, ncp2_dyad_mask, dna_ids, dimer, None, 'segname', ncp2_dyad_mol, debug)
            except RuntimeError:
                print 'ERROR: curve_fit failed for frame %d' % frame
                ncp2_origin = all_ncp2_origins[-1] *  0
                ncp2_axes = all_ncp2_axes[-1] * 0
        all_ncp2_origins.append(ncp2_origin)
        all_ncp2_axes.append(ncp2_axes)

    all_ncp2_axes = np.array(all_ncp2_axes)
    all_ncp2_origins = np.array(all_ncp2_origins)
    
    if aligned:
        all_ncp1_axes = np.copy(all_ncp2_axes)
        all_ncp1_origins = np.copy(all_ncp2_origins)
        all_ncp1_axes[:] = ncp1_axes
        all_ncp1_origins[:] = ncp1_origin
    else:
        all_ncp1_axes = np.array(all_ncp1_axes)
        all_ncp1_origins = np.array(all_ncp1_origins)
    
    return all_ncp1_origins, all_ncp2_origins, all_ncp1_axes, all_ncp2_axes

def print_min(magnitudes, mol, dist_to='origin'):
    i_min    = magnitudes.argmin()
    val_min  = magnitudes.min()
    coor_min = mol.coor()[0,i_min]
    print "min dist %s: %0.3f A, atom %d at [%0.3f, %0.3f, %0.3f]" % (
        dist_to, val_min, i_min, coor_min[0], coor_min[1], coor_min[2])
    print "%s: %s, resid %d, %s, %s" % (dist_to, mol.segname()[i_min], 
        mol.resid()[i_min], mol.resname()[i_min], mol.name()[i_min])
        
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
    
    # closest_atom_x_axis = [34710, 27887]
    # closest_atom_y_axis = [ 1563, 14817]
    # closest_atom_origin = [17061, 23178]
    
    main(dimer, dimer_dcdfile)

    
    print '\m/ >.< \m/'