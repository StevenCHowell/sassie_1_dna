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
    
def cylinder_distances(coor, R, X0, Y0, Vx, Vy):
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

def main(coor):
    ideal = np.zeros(len(coor))
    R = 42
    guess = (R, 0, 0, 1, 1)
    opt_params, cov_params= curve_fit(cylinder_distances, coor, ideal, guess)
    return opt_params

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
    R[:3,3] = origin
        
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

def show_cylinder(coor, params):
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
    array = np.linspace(-24,25)
    xs = array*Vx + X0
    ys = array*Vy + Y0
    zs = array*Vz + Z0
    ax.plot(xs, ys, zs, label='cylinder axis')
    
    ax.scatter(coor[:,0], coor[:,1], coor[:,2], label='data points')
    
    # X,Y,Z = cylinder(np.ones((10,1))*R, 20)
    X_raw, Y_raw, Z_raw = cylinder(R, 20)
    vector = np.array([Vx, Vy, Vz])
    origin = np.array([X0, Y0, Z0])
    X, Y, Z = transform_surf(X_raw, Y_raw, Z_raw, vector, origin)
    
    ax.plot_wireframe(X,Y,Z, label='cylinder')
    
    ax.legend()
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()
    
    print '\m/ >.< \m/'

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

if __name__ == '__main__':
    
    pdb_file = '1KX5tailfold_167bp.pdb'
    nuc = sasmol.SasMol(0)
    nuc.read_pdb(pdb_file)
    basis_filter = ' ( chain[i] ==  "I"  or chain[i] ==  "J"  ) and name[i] ==  "C1\'" '
    error, mask = nuc.get_subset_mask(basis_filter)
    error, coor = nuc.get_coor_using_mask(0, mask)
    params = main(coor[0]) # params = [R, X0, Y0, Vx, Vy]
    show_cylinder(coor[0], params)

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