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
from scipy.optimize import curve_fit
import logging
import sassie.sasmol.sasmol as sasmol

LOGGER = logging.getLogger(__name__) #add module name manually


class MainError(Exception):
    pass

def get_c1p_coors():
    NotImplemented
    
def cylinder_distances_from_R(coor, R, X0, Y0, Vx, Vy):
    origin = np.array([X0, Y0, 0]) # Z0 = 0

    Vz = 1
    length = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    Ux = Vx/length
    Uy = Vy/length
    Uz = Vz/length
    U = np.array([Ux,Uy,Uz])
    
    A = np.dot((coor-origin),U)              # component of array from origin to point along axis
    D = (coor-origin)-np.outer(A,U)          # vectors from axis to point
    dist = np.sqrt(np.square(D).sum(axis=1)) # distance from axis to point
    
    return dist-R

def vector_from_cylinder_axis(coor, R, X0, Y0, Vx, Vy):
    origin = np.array([X0, Y0, 0]) # Z0 = 0

    Vz = 1
    length = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    Ux = Vx/length
    Uy = Vy/length
    Uz = Vz/length
    U = np.array([Ux,Uy,Uz])
    
    coor_cyl = coor-origin
    A = np.dot((coor_cyl),U)     # component of array from origin to point along axis
    D = (coor_cyl)-np.outer(A,U) # vectors from axis to point
    
    return D

def transform_coor(coor3, vector, origin):
    # initialize vector arrays for coordinates and orientation vectors
    # changing them from 3 component vectors into 4 component vectors to
    # incorporate transaltions into the matrix math
    try:
        r, c = coor3.shape
        if c != 3:
            if r == 3:
                coor3 = coor3.T
    except:
        coor3 = coor3.reshape(1,3)
    
    coor4 = np.ones((len(coor3), 4))
    coor4[:, 0:3] = coor3     #; print 'coor4 =', coor4

    # angles to get to vector

    # create the translation-rotation matrix
    # This is intended to be multiplied from the right (unlike standard matrix
    # multiplication) so as not to require transposing the coordinate vectors.
    vector = vector.reshape(3)
    [vx, vy, vz] = vector

    tz = np.arctan2(vx, vy)
    cz = np.cos(tz)
    sz = np.sin(tz)

    tx = np.arctan2(vx*sz + vy*cz, vz)

    # initialize the rotation about X, then Z:
    Rx = rotate('x', -tx)
    Rz = rotate('z', -tz)
    Rxz = np.dot(Rx,Rz)
    R = np.eye(4)
    R[:3, :3] = Rxz
    R[3,:3] = origin
    
    
    # R would rotate V to [0, 0, 1], we want the reverse of that: V.T
    result = np.dot(coor4,R)

    return result[:,:-1]

def rotate(axis, theta):
    R = np.eye(3)
    ct = np.cos(theta)
    st = np.sin(theta)
    if axis.lower()=='x':
        (R[1, 1], R[1, 2]) = (ct, st)
        (R[2, 1], R[2, 2]) = (-st, ct)
    elif axis.lower()=='y':
        (R[0, 0], R[0, 2]) = (ct, -st)
        (R[2, 0], R[2, 2]) = (st, ct)
    elif axis.lower()=='z':
        (R[0, 0], R[0, 1]) = ( ct, st)
        (R[1, 0], R[1, 1]) = (-st, ct)
    else:
        assert True, "ERROR!!! did not recognize rotation axis"
    return R

def transform_surf(X, Y, Z, vector, origin):
    r, c = X.shape

    coor = np.array([X.reshape((r*c)),Y.reshape((r*c)),Z.reshape((r*c))]).T

    coor = transform_coor(coor, vector, origin)
    
    X_new = coor[:,0].reshape((r,c))
    Y_new = coor[:,1].reshape((r,c))
    Z_new = coor[:,2].reshape((r,c))

    return X_new, Y_new, Z_new

def show_cylinder(coor, params, nuc_origin, nuc_axes, dyad_origin):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    [R, X0, Y0, Vx, Vy] = params
    Z0 = 0
    Vz = 1
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # xs = [X0, X0+Vx]
    # ys = [Y0, Y0+Vy]
    # zs = [Z0, Z0+Vz]
    plot_cylinder_axis = False
    if plot_cylinder_axis:
        array = np.linspace(-24,25)
        xs = array*Vx + X0
        ys = array*Vy + Y0
        zs = array*Vz + Z0
        ax.plot(xs, ys, zs, label='cylinder axis')
    
    ax.plot(coor[:,0], coor[:,1], coor[:,2], 'o', label="C1' coordinates")

    plot_cylinder_origin = False
    if plot_cylinder_origin:
        # ax.scatter(X0, Y0, Z0, color='blue', label='origin', marker='x')
        cyl_origin = np.array([[X0, Y0, Z0]])
        # ax.plot(cyl_origin[:,0], cyl_origin[:,1], cyl_origin[:,2], color='blue', label='cylinder_origin', marker='x', line='')
        ax.plot(cyl_origin[:,0], cyl_origin[:,1], cyl_origin[:,2], 'bx', label='cylinder origin')
    
    # X,Y,Z = cylinder(np.ones((10,1))*R, 20)
    vector = np.array([Vx, Vy, Vz])
    X_raw, Y_raw, Z_raw = cylinder(np.ones((2,1))*R, 20)
    h = 40
    Z_raw = (Z_raw-0.5) * h
    X, Y, Z = transform_surf(X_raw, Y_raw, Z_raw, vector, nuc_origin)
    ax.plot_wireframe(X,Y,Z, label='cylinder', color='orange', lw='2')

    ax.plot(nuc_origin[:,0], nuc_origin[:,1], nuc_origin[:,2], 'gs', label='nucleosome origin') 
    styles = ['r-','g-','b-']
    labels = ['X-axis', 'Y-axis', 'Z-axis']
    for (i, axes) in enumerate(nuc_axes):
        axes_vec = np.concatenate((nuc_origin, axes*10 + nuc_origin), axis=0)
        ax.plot(axes_vec[:,0], axes_vec[:,1], axes_vec[:,2], styles[i], label=labels[i])

    dyad_origin = dyad_origin.reshape(1,3)
    ax.plot(dyad_origin[:,0], dyad_origin[:,1], dyad_origin[:,2], 'rs', label='dyad origin')

    # styles = ['r-','g-','b-']
    # labels = ['dyad X-axis', 'dyad Y-axis', 'dyad Z-axis']
    # for (i, axes) in enumerate(orientations):
        # dyad_axes = np.concatenate((dyad_origin, axes*10 + dyad_origin), axis=0)
        # ax.plot(dyad_axes[:,0], dyad_axes[:,1], dyad_axes[:,2], styles[i], label=labels[i])
    # dyad = dyad_mol.coor()[0]
    # ax.plot(dyad[:,0], dyad[:,1], dyad[:,2], 'ro', label='dyad bp')
    
    # ax.plot_wireframe(X_raw,Y_raw,Z_raw, label='cylinder', color='red')
    # ax.scatter(0, 0, 0, color='red', label='data points', marker='x')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.axis('equal')
    plt.legend(loc='upper left', numpoints=1)
    plt.show()
    
    print '\m/ >.< \m/'

def get_dna_bp_reference_frame(dna_ids, bp_mol, dna_id_type='segname'):
    '''
    The x-axis points in the direction of the major groove along what would 
    be the pseudo-dyad axis of an ideal Watson-Crick base-pair, i.e. the 
    perpendicular bisector of the C1'...C1' vector spanning the base-pair. 
    The y-axis runs along the long axis of the idealizosed base-pair in the
    direction of the sequence strand, parallel with the C1'...C1' vector, 
    and displaced so as to pass through the intersection on the 
    (pseudo-dyad) x-axis of the vector connecting the pyrimidine Y(C6) and 
    purine R(C8) atoms. The z-axis is defined by the right-handed rule, 
    i.e. z = x cross y. (doi:10.1006/jmbi.2001.4987)
    '''
    c6c8_string = ('(((resname[i] == "CYT" or resname[i] == "THY") and name[i] == "C6") or'
                   ' ((resname[i] == "GUA" or resname[i] == "ADE") and name[i] == "C8"))')

    dna1_c1p_filter  = '%s[i] == "%s" and name[i] == "C1\'" ' % (dna_id_type.lower(), dna_ids[0])
    dna2_c1p_filter  = '%s[i] == "%s" and name[i] == "C1\'" ' % (dna_id_type.lower(), dna_ids[1])
    dna1_c6c8_filter = '%s[i] == "%s" and %s' % (dna_id_type.lower(), dna_ids[0], c6c8_string)
    dna2_c6c8_filter = '%s[i] == "%s" and %s' % (dna_id_type.lower(), dna_ids[1], c6c8_string)

    e0, dna1_c1p_mask = bp_mol.get_subset_mask(dna1_c1p_filter)
    e1, dna2_c1p_mask = bp_mol.get_subset_mask(dna2_c1p_filter)
    e2, dna1_c6c8_mask = bp_mol.get_subset_mask(dna1_c6c8_filter)
    e3, dna2_c6c8_mask = bp_mol.get_subset_mask(dna2_c6c8_filter)
    assert np.sum(dna1_c1p_mask) == np.sum(dna1_c6c8_mask) == np.sum(dna2_c1p_mask) == np.sum(dna2_c6c8_mask), "ERROR: input did not contain atoms necessary for determining orientation"
    

    dna1_c1p = np.dot(dna1_c1p_mask, bp_mol.coor()[0])
    dna2_c1p = np.dot(dna2_c1p_mask, bp_mol.coor()[0])
    dna1_c6c8 = np.dot(dna1_c6c8_mask, bp_mol.coor()[0])
    dna2_c6c8 = np.dot(dna2_c6c8_mask, bp_mol.coor()[0])

    y_vec = dna1_c1p - dna2_c1p
    y_mag = np.sqrt(np.dot(y_vec, y_vec))
    y_hat = y_vec/y_mag
    
    # following: http://geomalgorithms.com/a05-_intersect-1.html
    Q0 = (dna1_c1p + dna2_c1p)/2
    P0 = dna2_c6c8
    P1 = dna1_c6c8
    w  = P0 - Q0
    u  = P1 - P0 
    s1 = np.dot(-y_hat, w)/np.dot(y_hat, u)
    assert 0 <= s1 <= 1, "ERROR: problem in calculating bead origin"
    bp_origin = P0 + s1 * u

    a = bp_origin - dna2_c1p
    x_vec = a - np.dot(a, y_hat) * y_hat
    x_hat = x_vec/np.sqrt(np.dot(x_vec, x_vec))
    
    z_hat = np.cross(x_hat, y_hat)
    
    bp_axes = np.array([x_hat, y_hat, z_hat])
    
    return bp_origin, bp_axes


def cylinder(r,n=20):
    '''
    Source: http://python4econ.blogspot.com/2013/03/matlabs-cylinder-command-in-python.html
    Returns the unit cylinder that corresponds to the curve r.
    INPUTS:  r - a vector of radii
             n - number of coordinates to return for each element in r

    OUTPUTS: x,y,z - coordinates of points
    '''

    # ensure that r is a column vector
    r = np.atleast_2d(r)
    r_rows,r_cols = r.shape

    # added to make it so the result is not just a circle
    if r_rows == r_cols == 1:
        r = np.ones((2,1))*r
    
    if r_cols > r_rows:
        r = r.T

    # find points along x and y axes
    points  = np.linspace(0,2*np.pi,n+1)
    x = np.cos(points)*r
    y = np.sin(points)*r

    # find points along z axis
    rpoints = np.atleast_2d(np.linspace(0,1,len(r)))
    z = np.ones((1,n+1))*rpoints.T
    
    return x,y,z

def get_dna_bp_and_axes(dna_resids, dna_ids, dna_mol):
    bp_filter = '( chain[i] == "%s" and resid[i] == %d ) or ( chain[i] == "%s" and resid[i] == %d )' % (dna_ids[0], dna_resids[0], dna_ids[1], dna_resids[1])
    error, bp_mask = dna_mol.get_subset_mask(bp_filter)
    bp_mol = sasmol.SasMol(0)
    error = dna_mol.copy_molecule_using_mask(bp_mol, bp_mask, 0)
    bp_origin, bp_axes = get_dna_bp_reference_frame(dna_ids, bp_mol, 'chain')

    return bp_origin, bp_axes, bp_mol

def get_ncp_origin_and_axes(ncp_c1p_filter, dyad_dna_resids, dyad_dna_id, ncp, plt_res=False):
    error, c1p_mask = ncp.get_subset_mask(ncp_c1p_filter)
    error, coor = ncp.get_coor_using_mask(0, c1p_mask)
    coor = coor[0]
    ideal = np.zeros(len(coor))
    R = 42
    guess = (R, 0, 0, 1, 1)
    
    ## fit a cylinder
    import time
    tic = time.time()
    opt_params, cov_params= curve_fit(cylinder_distances_from_R, coor, ideal, guess)
    toc = time.time() - tic
    print 'fitting a cylinder to the NCP took %0.3f s' %toc
    
    [R, X0, Y0, Vx, Vy] = opt_params
    Z0 = 0
    Vz = 1
    cyl_origin = np.array([X0, Y0, Z0])
    z = np.array([Vx, Vy, Vz])
    z_hat = z/np.sqrt(np.dot(z,z))
    
    dyad_origin, dyad_axes, dyad_mol = get_dna_bp_and_axes(dyad_dna_resids, dyad_dna_id, ncp)

    ## calculate distance from dyad_orign to the axis
    x = vector_from_cylinder_axis(dyad_origin, opt_params[0], opt_params[1], opt_params[2], opt_params[3], opt_params[4])
    ncp_origin = dyad_origin - x

    # xp0 = vector_from_cylinder_axis(dyad_origin, params[0], params[1], params[2], params[3], params[4])
    # xp1 = xp0-origin
    # x = xp1 - np.dot(xp1, z_hat)*z_hat #subtract from x the projection along z_hat    
    
    x = x.reshape(3)
    x_hat = x/np.sqrt(np.dot(x,x))
    
    y_hat = np.cross(z_hat,x_hat)
    ncp_axes = np.array([x_hat, y_hat, z_hat])

    if plt_res:
        ## display the fit results
        show_cylinder(coor, opt_params, ncp_origin, ncp_axes, dyad_origin)
    
    return ncp_origin, ncp_axes

if __name__ == '__main__':
    import time
    pdb_file = '1KX5tailfold_167bp.pdb'
    ncp = sasmol.SasMol(0)
    ncp.read_pdb(pdb_file)
    basis_filter = ' ( chain[i] ==  "I"  or chain[i] ==  "J"  ) and name[i] ==  "C1\'" '

    dyad_dna_resids = [0, 0]
    dyad_dna_id = ['I', 'J']
    tic = time.time()
    ncp_origin, ncp_axes = get_ncp_origin_and_axes(basis_filter, dyad_dna_resids, dyad_dna_id, ncp) 
    toc = time.time() - tic
    print 'determining the NCP origin and axes took %0.3f s' %toc
    print 'ncp_origin =', ncp_origin
    print 'ncp_axes =\n', ncp_axes

    # coor = np.array([[0,0,1]])    
    # vector = np.array([params[3], params[4], 1])
    # origin = np.array([params[2], params[3], 0])
    # # v = np.ones(3)
    # # origin = np.zeros(3)
    # res = transform_coor(coor, vector, origin)
    # print '|res|:', res/np.sqrt(np.dot(res,res))
    # print 'should parrallel'
    # print '|v|:', vector/np.sqrt(np.dot(vector,vector))
    
    print '\m/ >.< \m/'