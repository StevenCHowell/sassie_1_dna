#!/usr/bin/env python
#coding:utf-8
#
# Author:  Steven C. Howell
# Purpose: Store the methods removed from the cgDNA_move.py module
# Created: 28 August 2014
# $Id$
#

import sys
import os
import os.path as op
import subprocess
import logging
import argparse
import numpy

LOGGER = logging.getLogger(__name__) #add module name manually


class MainError(Exception):
    pass


def write_xyz(filename, coor, frame):

    natoms=len(coor)
    if(frame==0):
        outfile=open(filename, 'w')
    else:
        outfile=open(filename, 'a')

    #outfile.write("%i\n" % natoms)
    #outfile.write('%s\t%s\n' % ('B', str(comment)))
    for i in range(natoms):
        #outfile.write('%s\t%f\t%f\t%f\n' % ('C', coor[i][0], coor[i][1],
        #                                    coor[i][2]))
        outfile.write('%f\t%f\t%f\n' % (coor[i][0], coor[i][1], coor[i][2]))
    outfile.close()

def make_bead_model(all_atom_pdb):
    '''
    this function creates a bead model from DNA
    this is the original function created by JC
    make_model will replace this by taking specific DNA to coarse grain as input
    '''
    aa_dna = sasmol.SasMol(0)
    aa_dna.read_pdb(all_atom_pdb)
    natoms = aa_dna.natoms()

    # print 'natoms = ', natoms

    cg_dna = sasmol.SasMol(0)
    basis_filter = ("(resname[i] == 'ADE' or resname[i] == 'GUA') and "
                    "(name[i] == 'N1' and (resid[i] < 21) "
                    "and segname[i]=='DNA1')")
    error, mask = aa_dna.get_subset_mask(basis_filter)

    frame = 0
    error = aa_dna.copy_molecule_using_mask(cg_dna, mask, frame)
    error, coor=aa_dna.get_coor_using_mask(frame, mask)

    cg_dna.setCoor(coor)

    diameter = 6.0*3.4

    cg_natoms = cg_dna.natoms()
    #        print 'cg_natoms = ', cg_natoms
    new_coor = numpy.zeros((1, cg_natoms, 3), numpy.float)
    #        print new_coor

    for i in xrange(cg_natoms):
        #                print "i = ", i
        new_coor[0][i][2] = i*diameter

    cg_dna.setCoor(new_coor)

    index = cg_dna.index()

    infile=open('dum.pdb', 'w')
    for i in xrange(cg_dna.natoms()):
        this_index = index[i]

        sx = cg_dna._coor[frame, i, 0]#[:8]
        sy = cg_dna._coor[frame, i, 1]#[:8]
        sz = cg_dna._coor[frame, i, 2]#[:8]

        infile.write(("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      "
            "%-4s%2s%2s\n") % (cg_dna.atom()[i], this_index, cg_dna.name()[i],
            cg_dna.loc()[i], cg_dna.resname()[i], cg_dna.chain()[i],
            cg_dna.resid()[i], cg_dna.rescode()[i], sx, sy, sz,
            cg_dna.occupancy()[i], cg_dna.beta()[i], cg_dna.segname()[i],
            cg_dna.element()[i], cg_dna.charge()[i]))

    infile.write('END\n')
    infile.close()
    os.remove('dum.pdb') # remove the temporary pdb file

    cg_dna.write_pdb("cg_test.pdb", frame, 'w')

    return cg_dna

def checkMag(vec):
    '''
    module to check that the magnitude of all vectors is the same
    '''
    erLim = 1e-3
    (r, c) = vec.shape
    m = numpy.zeros(r)
    for i in xrange(r):
        m[i] = mag(vec[i])

    avMag = numpy.mean(m)
    test = numpy.abs(m-avMag) > erLim  # check for differences in l-lengths
    if test.any():
        print 'ERROR: the vectors do not have uniformly magnitude'
        print 'm = \n', m
        print test
        print 'vec = ', vec

    return (avMag)

def p_energy_wca_mixed(w_d, coor_d, wca0_d, trial_bead, w_p, p_coor,
                       wca0_p, group_mask_p, wca0_dp):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    '''
    w_dp = numpy.mean([w_d, w_p])

    wca1_d = numpy.copy(wca0_d)
    wca1_dp = numpy.copy(wca0_dp)
    wca1_p = numpy.copy(wca0_p)

    test_d = 2.**(1./6.)*w_d
    test_p = 2.**(1./6.)*w_p
    test_dp = 2.**(1./6.)*w_dp

    (N_d, col) = d_coor.shape
    (N_p, col) = p_coor.shape

    # calculate the distance between moved DNA beads and all other DNA beads
    for i in xrange(trial_bead, N_d):    # only calculate changed distances
        ri = d_coor[i, :]                # get ith bead's coordinates
        for j in xrange(i):
            rj = d_coor[j, :]            # get the coordinates for the jth bead
            # calculate the distance between the ith and jth beads
            rij = numpy.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            if rij < test_d:
                wca1_d[i, j] = (w_d/rij)**12.-(w_d/rij)**6.+0.25

    # calculate the distance between DNA beads and proteins
    for i in xrange(N_d):
        ri = d_coor[i, :]
        for j in xrange(N_p):
            rj = p_coor[j, :]
            # get the distance between the ith DNA and jth protein beads
            rij = numpy.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            if rij < test_dp:
                wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25

    # calculate the distance between proteins and proteins
    for i in xrange(N_p):
        ri = p_coor[i, :]
        for j in xrange(i):
            rj = p_coor[j, :]
            # get the distance between the ith DNA and jth protein beads
            rij = numpy.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            if rij < test_p:
                wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25

    res_d = 4 * numpy.sum(wca1_d)
    res_p = 4 * numpy.sum(wca1_p)
    res_dp = 4 * numpy.sum(wca1_dp)

    res = res_d + res_p + res_dp
    return (res, wca1_d, wca1_p, wca1_dp)

def p_energy_wca_mixed_fast(w_d, coor_d, wca0_d, trial_bead, w_p, p_coor,
                             wca0_p, group_mask_p, wca0_dp, initial=False):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads

    this method skips unnecessary distances, (compared to the longer method,
    the only difference in result could be attributed to float error)
    '''
    w_dp = numpy.mean([w_d, w_p])

    wca1_d = numpy.copy(wca0_d)
    wca1_dp = numpy.copy(wca0_dp)
    wca1_p = numpy.copy(wca0_p)

    test_d = 2.**(1./6.)*w_d
    test_p = 2.**(1./6.)*w_p
    test_dp = 2.**(1./6.)*w_dp

    (N_d, col) = d_coor.shape
    (N_p, col) = p_coor.shape

    # calculate the distance between moved DNA beads and all other DNA beads
    for i in xrange(trial_bead, N_d):    # only calculate changed distances
        ri = d_coor[i, :]                # get ith bead's coordinates
        for j in xrange(i):
            rj = d_coor[j, :]            # get the coordinates for the jth bead
            # calculate the distance between the ith and jth beads
            rij = numpy.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            if rij < test_d:
                wca1_d[i, j] = (w_d/rij)**12.-(w_d/rij)**6.+0.25

    # get the indices of the moved and stationary proteins
    ind_moved_p = mask2ind(group_mask_p)
    ind_stationary_p = mask2ind(-(group_mask_p-1))

    if initial:
        for i in xrange(N_d):           # all the DNA
            ri = d_coor[i, :]           # get ith bead's coordinates
            for j in xrange(N_p):       # all the protein
                rj = p_coor[j, :]
                # get the distance between the ith DNA and jth protein beads
                rij = numpy.sqrt((ri[0]-rj[0])**2. +
                              (ri[1]-rj[1])**2. +
                              (ri[2]-rj[2])**2.)
                if rij < test_dp:
                    wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25
        for i in xrange(N_p):
            ri = p_coor[i, :]
            for j in xrange(i):
                rj = p_coor[j, :]
                # get the distance between the ith DNA and jth protein beads
                rij = numpy.sqrt((ri[0]-rj[0])**2. +
                              (ri[1]-rj[1])**2. +
                              (ri[2]-rj[2])**2.)
                if rij < test_p:
                    wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25
    else:
        # calculate the distance between moved DNA and stationary proteins
        if len(ind_stationary_p) > 0:
            for i in xrange(trial_bead, N_d): # only the DNA that moved
                ri = d_coor[i, :]             # get ith bead's coordinates
                for j in ind_stationary_p:    # only the stationary proteins
                    rj = p_coor[j, :]
                    # get distance between the ith DNA and jth protein beads
                    rij = numpy.sqrt((ri[0]-rj[0])**2. +
                                  (ri[1]-rj[1])**2. +
                                  (ri[2]-rj[2])**2.)
                    if rij < test_dp:
                        wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25
        # calculate the distance between stationary DNA and moved proteins
        if len(ind_moved_p) > 0:
            for i in xrange(trial_bead): # only want the DNA that did not move
                ri = d_coor[i, :]        # get ith bead's coordinates
                for j in ind_moved_p:    # only want stationary proteins
                    rj = p_coor[j, :]
                    # get distance between the ith DNA and jth protein beads
                    rij = numpy.sqrt((ri[0]-rj[0])**2. +
                                  (ri[1]-rj[1])**2. +
                                  (ri[2]-rj[2])**2.)
                    if rij < test_dp:
                        wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25

        # calculate the distance between moved proteins and stationary proteins
        if len(ind_moved_p) > 0 and len(ind_stationary_p) > 0:
            for i in ind_moved_p:
                ri = p_coor[i, :]
                for j in ind_stationary_p:
                    rj = p_coor[j, :]
                    # get distance between the ith DNA and jth protein beads
                    rij = numpy.sqrt((ri[0]-rj[0])**2. +
                                  (ri[1]-rj[1])**2. +
                                  (ri[2]-rj[2])**2.)
                    if rij < test_p:
                        wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25

    res_p = 4 * numpy.sum(wca1_p)
    res_d = 4 * numpy.sum(wca1_d)
    res_dp = 4 * numpy.sum(wca1_dp)

    res = res_d + res_p + res_dp
    return (res, wca1_d, wca1_p, wca1_dp)

def f_energy_wca_mixed(w_d, coor_d, wca0_d, trial_bead0, w_p, p_coor, wca0_p,
                     group_mask_p, wca0_dp, initial=False):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN

    this method skips unnecessary distances, compared to other method,
    difference in results could be attributed to float error
    '''
    #import sys ;
    #sys.path.append('./')
    #sys.path.append('/home/schowell/Dropbox/gw_phd/code/pylib/sassie/')

    w_dp = numpy.mean([w_d, w_p])
    wca1_d = numpy.copy(wca0_d)
    wca1_dp = numpy.copy(wca0_dp)
    wca1_p = numpy.copy(wca0_p)

    (N_d, col) = d_coor.shape
    (N_p, col) = p_coor.shape

    # get the indices of the moved and stationary proteins
    ind_moved_p = mask2ind(group_mask_p)
    ind_stationary_p = mask2ind(-(group_mask_p-1))
    if len(ind_moved_p):
        imoved0 = numpy.min(ind_moved_p)
    else:
        imoved0 = len(ind_stationary_p) - 1
    imoved1 = imoved0 + 1
    trial_bead1 = trial_bead0 + 1

    # calculate the distance between moved DNA beads and all other DNA beads
    wca1_d = collision.wca_d(d_coor, trial_bead1, w_d, wca1_d)
    # trial bead index will begin at 1 instead of 0

    res_d = 4.*numpy.sum(wca1_d)

    if initial:
        wca1_p = collision.wca_d(p_coor, 1, w_p, wca1_p)
        res_p = 4 * numpy.sum(wca1_p)

        wca1_dp = collision.wca_nbym(d_coor, p_coor, 1, N_d, 1, N_p, w_dp,
                                     wca1_dp)
        res_dp = 4 * numpy.sum(wca1_dp)

    else:
        # calculate the distance between moved proteins and stationary proteins
        if len(ind_moved_p) > 0 and len(ind_stationary_p) > 0:
            wca1_p = collision.wca_nbyn(p_coor, imoved1, N_p, 1, imoved0, w_p,
                        wca1_p) # trial bead index will begin at 1 instead of 0
        res_p = 4 * numpy.sum(wca1_p)

        # calculate the distance between moved DNA and stationary proteins
        if len(ind_stationary_p) > 0:
            wca1_dp = collision.wca_nbym(d_coor, p_coor, trial_bead1, N_d, 1,
                        imoved0, w_dp, wca1_dp)
        # calculate the distance between stationary DNA and moved proteins
        if len(ind_moved_p) > 0:
            wca1_dp = collision.wca_nbym(d_coor, p_coor, 0, trial_bead0,
                        imoved1, N_p, w_dp, wca1_dp)
        res_dp = 4 * numpy.sum(wca1_dp)

    # FIX THIS:
    #s res = res_d + res_p + res_dp # res_dp is too high

    res = res_d + res_p # leave out DNA-protein interaction energy
    return (res, wca1_d, wca1_p, wca1_dp)

def p_energy_wca(w, coor, wca0, trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    '''
    wca1 = numpy.copy(wca0)
    test = 2.**(1./6.)*w
    (N, col) = coor.shape
    for i in xrange(trial_bead, N):
        ri = coor[i, :]
        for j in xrange(0, i):
            rj = coor[j, :]
            rij = numpy.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            #s print '(2^(1/6)*w, rij) = ', (test, rij)
            if rij < test:
                wca1[i, j] = (w/rij)**12.-(w/rij)**6.+0.25
    res = 4*numpy.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def test_wca(val1, val2, out=False):
    if val1.shape != val2.shape:
        print 'ERROR: val1 and val2 are not comparable'
        return

    diff = (val1 - val2)
    per_diff = diff / val1
    if numpy.any(numpy.isnan(per_diff)):
        per_diff[numpy.isnan(per_diff)] = 0
    sum_diff = numpy.sum(per_diff)
    max_diff = numpy.max(numpy.abs(per_diff))

    print 'val1:\n', val1
    print 'val2:\n', val2
    print 'difference:\n', diff
    print 'percent difference:\n', per_diff
    print 'total percent difference:\n', sum_diff
    print 'maximum percent difference:\n', max_diff

    if out:
        return (per_diff, max_diff, diff)

    return ''

def mag(vec):

    (r, ) = vec.shape
    sumSquares = 0.0
    for i in xrange(r):
        sumSquares += vec[i]**2

    return numpy.sqrt(sumSquares)

def bead2bp_com(d_coor, vecZ, bps, dna_type, use_ru, bp_coor):
    '''
    Given a bead that represents several bps, create a bead for each bp
    Was going to use this for the collision calculator but it was not necessary

    Usage:
    bp_coor = numpy.zeros((d_nbeads*3, 3))
    bead2bp_com(d_coor, vecXYZ[2], ARGS.bp_per_bead, 'b', True, bp_coor)
    print bp_coor

    '''
    if dna_type.lower() == 'b':
        s = 3.4
    elif dna_type.lower() == 'a':
        s = 2.56
    elif  dna_type.lower() == 'z':
        s = 3.7
    else:
        print '\n>>> ERROR: bad DNA type given <<<\n'
        s = 0

    assert len(d_coor) == len(vecZ), 'mismatched N coordinates and Z-vectors'
    for bead in xrange(len(d_coor)): # only the DNA that moved
        dna_i = d_coor[bead, :]          # get bead's coordinates
        for bp in xrange(bps):
            offset = s * ((bps - 1.0) / 2.0 - bp)
            if bead > len(u) - 1:
                z = u[bead-1, :]
            elif offset < 0 and bead > 0:
                z = u[bead-1, :]
            else:
                z = u[bead, :]

            z = z / numpy.sqrt(z[0]**2+z[1]**2+z[2]**2)
            if use_ru:
                ru = dna_i + z * offset
                bp_coor[bead * bps + bp] = ru
            else:
                rz = dna_i + vecZ[bead] * offset
                bp_coor[bead * bps + bp] = rz

    '''
    # Code to generate plot of coordinates:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    ax = Axes3D(plt.figure())
    ax.scatter(d_coor[:,0], d_coor[:,1], d_coor[:,2],color='red',s=50)
    ax.scatter(bp_coor[:,0], bp_coor[:,1], bp_coor[:,2],color='green')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    r_max = 15
    ax.set_xlim([-r_max,r_max])
    ax.set_ylim([-r_max,r_max])
    plt.show()
    '''
def get_distances(cg_dna, cg_pro, pdb):
    write_xyz(ARGS.pdb[:-4]+'_cgDNA.dat', cg_dna.coor()[0], 0)
    write_xyz(ARGS.pdb[:-4]+'_cgPro.dat', cg_pro.coor()[0], 0)

    import survey_distances as sd
    dna_file_name = pdb[:-4] + '_cgDNA.pdb'
    pro_file_name = pdb[:-4] + '_cgPro.pdb'
    cg_dna.write_pdb(dna_file_name, 0, 'w')
    cg_pro.write_pdb(pro_file_name, 0, 'w')
    dna_pro_dist = sd.f_calc_dist2(cg_dna.coor()[0], cg_pro.coor()[0])
    numpy.savetxt('dna-pro-dist_c11_withTails.txt', dna_pro_dist)
    mn = min(dna_pro_dist)
    mx = max(dna_pro_dist)
    print 'min distance ', mn
    print 'max distance ', mx
    sd.plot_historam(dna_pro_dist,'DNA-Protein Distances', int(1.1*mx))

def centerPDB(pdb_file):
    '''
    load and center given pdb file
    '''
    frame = 0
    pdb = sasmol.SasMol(frame)
    pdb.read_pdb(pdb_file)    
    pdb.center(frame)
    out_pdb_name = pdb_file[0:-4] + "_centered.pdb"
    print 'created centered pdb: ', out_pdb_name
    pdb.write_pdb(out_pdb_name, frame, 'w')

def makeLongDNA(n_lp):
    print 'making DNA that is %d*lp long' %n_lp

    # 15 bp/bead or 51 A/bead (3.4 A/bp)

    lp = 530 # persistence length in A
    l = 2**(1./6.)*46   # separation distance between beads = 51.6A

    longDNA = sasmol.SasMol(0)
    L = n_lp*lp
    N = int(L/l)
    natoms = N+1
    print 'natoms = ', natoms
    longDNA._L = L
    longDNA._natoms = natoms

    longCoor = numpy.zeros((1, natoms, 3), numpy.float) # initialize long DNA position
    longCoor[0][:, 2] = range(natoms) # set the z-values to the array index
    longCoor *= l                   # scale the z-values to the right seperation
    # print longCoor[-5:]

    longDNA.setCoor(longCoor)
    longDNA.setElement(['C']*natoms)

    vecXYZ = numpy.zeros((3, natoms, 3))
    vecXYZ[0] = [1, 0, 0]
    vecXYZ[1] = [0, 1, 0]
    vecXYZ[2] = [0, 0, 1]
    # n = L/l                         # number of times need to repeat the grain
    # print '(l, L, n)', (l, L, n)

    return (longDNA, vecXYZ)

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'Generate modified DNA or DNA-protein structures'
        #epilog = 'no epilog found'
    )
    parser.add_argument("-s", "--seed", default=0, type=int,
        help = ("Seed for generating random number for go-back routine "
                "(default: use system time)") )
    parser.add_argument("-gb", "--goback", default=-1, type=int,
        help = ("Number of fails before going back to previously accepted "
                "structure. Input should be > 0, default = -1 (off)") )
    parser.add_argument("-nd", "--n_dcd_write", default=1, type=int,
        help = ("Number of Monte Carlo steps before saving a dcd step; "
              "default = 1"))
    parser.add_argument("-n", "--nsteps", default=100, type=int,
        help = ("Number of Monte Carlo steps to perform; default = 100"))
    #    parser.add_argument("-s", "--show", action="store_true",
    #        help="show a plot of the configuration after each iteration")
    parser.add_argument("-tm", "--theta_max", nargs='+', type=numpy.float,
        help = ("Max rotation angle for the x and y rotation; "
                "max z scales with this paramater (curently 1x but code is "
                "designed to allow max z to scale differently from max x and "
                "y); default = 5"), default=[25])
    parser.add_argument("-ns", "--n_soft", default=1, type=int,
        help = ("Number of bead over which to apply any rotation (i.e. soften "
                "the bending); default = 1") )
    parser.add_argument("-bp", "--bp_per_bead", default=1, type=int,
        help = "Number of bp for each coarse-grained bead; default = 1")
    parser.add_argument("-t", "--temperature", default=300, type=float,
        help = "Temperature for calculating the coulomb energy; default = 278")
    parser.add_argument("-fl", "--flex", default=False, action="store_true",
        help = ("Generate output file with the flexible resids formatted into "
                "2 columns.  The first row contains the segnames for the 2 "
                "DNA strands; default = False"))

    parser.add_argument("-rp", "--rm_pkl", default=False, action="store_true",
                        help="remove the pkl file containing cg parameters")

    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="display the acceptance rate and debug info")

    parser.add_argument("-kcg", "--keep_cg_files", default=False, 
                        action="store_true",
                        help="keep the coarse-grained output files")

    parser.add_argument("-u", dest="keep_unique", action="store_true",
        help="only store coordinates of unique structures; the default")
    parser.add_argument("-a", dest="keep_unique", action="store_false",
        help="store coordinates for all accepted structures")
    parser.set_defaults(keep_unique=True)

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-l", "--Llp", type=int,
        help=("length of desired coarse-grain long DNA to simulate (in units of"
              " persistence lengths); -->THIS FUNCTIONALITY NO LONGER WORKS "
              "(NEED TO REVISE FOR BACKWARD COMPATIBILITY)<--") )
    group2.add_argument("-p", "--pdb", default="new_dsDNA.pdb",
        help="all atom pdb file, default='new_dsDNA.pdb'")

    return parser.parse_args()

def main():
    NotImplemented    

if __name__ == '__main__':
    import unittest    
    
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    args = parse()  # this makes args global
    # if args.pdb is not None:
        # args.pdb = os.getcwd() + '/' + args.pdb

    s_theta_max = ''
    for (i, theta) in enumerate(args.theta_max):
        if i > 0:
            s_theta_max += ', '
        s_theta_max += str(theta)


    main()
    unittest.main()
