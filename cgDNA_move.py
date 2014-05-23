#!/usr/bin/evn python
# Author:   --<Steven Howell>
# Purpose:  Provide structure movement using SASSIE protocols
# Created: 12/01/2013
# $Id: cgDNA_move.py,v 1.27 2014-05-23 19:42:20 schowell Exp $

# time using FORTRAN double loop, N=1000, iters=1000 (so 1*10^6 steps): 958.887075186 seconds
# time using python double loop, N=1000, iters=1000 (so 1*10^6 steps):

import sassie.sasmol.sasmol as sasmol
import numpy as np
import string
import os
import locale
import sys
import random
import os.path as op
import subprocess
import logging
import group_rotation
import sys ; sys.path.append('./')
#import sys ; sys.path.append('/home/schowell/Dropbox/gw_phd/code/pylib/sassie/')
import collision
import time

LOGGER = logging.getLogger(__name__) #add module name manually

class MainError(Exception):
    pass

def write_xyz(filename, coor, comment, frame):

    natoms=len(coor)
    if(frame==0):
        outfile=open(filename, 'w')
    else:
        outfile=open(filename, 'a')

    outfile.write("%i\n" % natoms)
    outfile.write('%s\t%s\n' % ('B', str(comment)))
    for i in range(natoms):
        outfile.write('%s\t%f\t%f\t%f\n' % ('C', coor[i][0], coor[i][1], coor[i][2]))
    outfile.close()

    return


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
    basis_filter = "(resname[i] == 'ADE' or resname[i] == 'GUA') and (name[i] == 'N1' and (resid[i] < 21) and segname[i]=='DNA1')"
    error, mask = aa_dna.get_subset_mask(basis_filter)

    frame = 0
    error = aa_dna.copy_molecule_using_mask(cg_dna, mask, frame)
    error, coor=aa_dna.get_coor_using_mask(frame, mask)

    cg_dna.setCoor(coor)

    diameter = 6.0*3.4

    cg_natoms = cg_dna.natoms()
    #        print 'cg_natoms = ', cg_natoms
    new_coor = np.zeros((1, cg_natoms, 3), np.float)
    #        print new_coor

    for i in xrange(cg_natoms):
        #                print "i = ", i
        new_coor[0][i][2] = i*diameter

    cg_dna.setCoor(new_coor)
    #        print new_coor
    #        print cg_dna.coor()

    #        charge = cg_dna.beta()
    #        cg_dna.setCharge(charge)

    #        print 'loc = ', cg_dna.loc()
    #        print 'rescode = ', cg_dna.rescode()
    #        print 'charge = ', cg_dna.charge()

    #        print 'type(index) = ', type(cg_dna.index())
    index = cg_dna.index()

    #sh        coor = cg_dna.coor()[0]
    #sh        comment = 'test'
    #sh        write_xyz('cg_dna.xyz', coor, comment, frame)

    infile=open('dum.pdb', 'w')
    for i in xrange(cg_dna.natoms()):
        this_index = index[i]

        sx = cg_dna._coor[frame, i, 0]#[:8]
        sy = cg_dna._coor[frame, i, 1]#[:8]
        sz = cg_dna._coor[frame, i, 2]#[:8]

        infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s\n" % (cg_dna.atom()[i], this_index, cg_dna.name()[i], cg_dna.loc()[i], cg_dna.resname()[i], cg_dna.chain()[i], cg_dna.resid()[i], cg_dna.rescode()[i], sx, sy, sz, cg_dna.occupancy()[i], cg_dna.beta()[i], cg_dna.segname()[i], cg_dna.element()[i], cg_dna.charge()[i]))

    infile.write('END\n')
    infile.close()
    os.remove('dum.pdb') # remove the temporary pdb file

    cg_dna.write_pdb("cg_test.pdb", frame, 'w')

    #sh        cg_dna.write_dcd("cg_test.dcd")

    return cg_dna

def make_cg_model(all_atom_pdb, dna_chains, dna_resids, l, pro_groups,
                  bp_perBead):
    '''
    this function creates a cg bead model from DNA chains supplied as input
    make_cgDNA_model takes specific DNA to coarse grain as input
    '''
    # load in the all atom pdb
    aa_all = sasmol.SasMol(0)
    aa_all.read_pdb(all_atom_pdb)
    natoms = aa_all.natoms() #; print 'natoms = ', natoms
    print 'total atoms = ', natoms
    frame = 0

    #~~~ DNA ONLY SECTION ~~~#
    chain1 = dna_chains[0]
    chain2 = dna_chains[1]
    resid1 = dna_resids[0]
    resid2 = dna_resids[1]
    # check the input
    assert np.abs(resid1[1]-resid1[0]) == np.abs(resid2[1]-resid2[0]), "number of paired bases in chains are not the same"

    bps = np.abs(resid1[1]-resid1[0])+1
    nbeads = int(np.round(bps/float(bp_perBead)))
    #print 'nbeads = ', nbeads

    # load the dna into its own sasmol object
    aa_dna = sasmol.SasMol(0)
    error, mask = aa_all.get_subset_mask('((moltype[i] == "dna" or moltype[i] == "rna"))')  # THIS SHOLUD NOT HAVE or rna
    error = aa_all.copy_molecule_using_mask(aa_dna, mask, frame) # pick out the dna
    aa_dna.write_pdb("aaDNA.pdb", frame, 'w')
    natomsD = aa_dna.natoms()
    print 'dna atoms = ', natomsD

    # populate the cg_dna sasmol object with atom properties
    cg_dna = sasmol.SasMol(0)
    basis_filter = "((name[i] == 'N1') and (resid[i] <= "+str(nbeads)+") and (chain[i] == '"+chain1+"'))" #; print 'bf = ', basis_filter
    error, mask = aa_dna.get_subset_mask(basis_filter)   #; print mask
    s = np.sum(mask)
    assert s == nbeads, "\n>>> ERROR!!! wrong number of beads created: expected %d but found %d\n" %(nbeads, s)

    error = aa_dna.copy_molecule_using_mask(cg_dna, mask, frame) # pick out the first N=nbeads 'P' atoms
    cg_coor = np.zeros((1, nbeads, 3)) # initialize array to store the coordinates of the beads

    r1a = resid1[0]
    r2a = resid2[0]
    all_beads = []    # list to contain all the all atom sasmol object of each bead
    dna_masks = []       # list to contain all the all atom masks
    link = np.zeros((nbeads-1, 2), dtype=int)  # these respresent the links between beads, except on the end, if a bead is flexible both it's links are flexible

    # to make flexid1 an array of the indices of flexible residues
    flexResids = [item for sublist in l for item in sublist]
    groupid = np.zeros(len(flexResids))
    # groupid = np.zeros(sum(len(x) for x in l))
    n_start = 0
    for i in xrange(len(l)):
        n_end = len(l[i]) + n_start
        groupid[n_start:n_end] = np.ones(len(l[i])) * i
        n_start = n_end

    flexResids = np.concatenate(([flexResids], [groupid]),0).T

    print "coarse-graining the DNA, this may take awhile..."
    for j in xrange(nbeads):
        if j+1 == nbeads:
            bp_perBead = bps - bp_perBead*j  # to account for the remainder

        bead = sasmol.SasMol(0)

        # Get the atoms from DNA chain 1
        if resid1[0] < resid1[1]:
            r1b = r1a + bp_perBead  # r1b > r1a
            basis_filter1 = "((resid[i] > "+str(r1a-1)+" and resid[i] < "+str(r1b)+") and (chain[i]=='"+chain1+"')) or "
        else:
            r1b = r1a - bp_perBead  # r1b < r1a
            basis_filter1 = "((resid[i] > "+str(r1b)+" and resid[i] < "+str(r1a+1)+") and (chain[i]=='"+chain1+"')) or "

        # Get the atoms from DNA chain 2
        if resid2[0] < resid2[1]:
            r2b = r2a + bp_perBead  # r2b > r2a
            basis_filter2 = "((resid[i] > "+str(r2a-1)+" and resid[i] < "+str(r2b)+") and (chain[i]=='"+chain2+"'))"
        else:
            r2b = r2a - bp_perBead   # r2b < r2a
            basis_filter2 = "((resid[i] > "+str(r2b)+" and resid[i] < "+str(r2a+1)+") and (chain[i]=='"+chain2+"'))"

        #s print '(r1a, r1b, r2a, r2b) = ' , (r1a, r1b, r2a, r2b)
        basis_filter = basis_filter1+basis_filter2     #; print 'basis_filter = ', basis_filter

        error, mask = aa_dna.get_subset_mask(basis_filter)   # create a mask to select the atoms for the bead
        dna_masks.append(mask)   # store the mask for the reverse coarse-graining
        error = aa_dna.copy_molecule_using_mask(bead, mask, frame)

        #s name = 'bead'+str(j)+'.pdb' ; bead.write_pdb(name, 0, 'w')                 # output each bead to a pdb for verification purposes

        com = bead.calccom(frame)
        cg_coor[0, j, :] = com  #; print 'com ', j, '= ', com       # store the coordinates of the com for the bead coordinate

        bead.center(0)  #; print   'bead._com = ', bead._com    # calculate the coodinates of the atoms in the bead with reference to the bead com
        all_beads.append(bead)

        # check residues from chain 1 to see if this is a flexible bead
        for i in xrange(r1a, r1b):
            if i in flexResids[:, 0]:
                (r, c) = np.where(flexResids==i)
                # this stores the assigned group in the temp var: beadGroup
                beadGroup = flexResids[r[0], 1]
                if j==0:                         # first bead
                    link[j] = [1, beadGroup]
                elif j==nbeads-1:                # last bead
                    link[j-2] = [1, beadGroup]
                else:  #all other beads, this is only 2 links in python
                    link[j-1:j+1] = [1, beadGroup]
                break  # bead defined as flexible, stop looping over residues

        # setup for next iteration
        r1a = r1b
        r2a = r2b

    # set the bead coordinates to the calculated com coordinates
    #s print 'cg_coor = \n', cg_coor
    beadgroups = np.zeros(len(link[:, 1])+1, dtype=int)
    beadgroups[1:] = link[:, 1]
    # need a trial bead for every flexible link (0th bead will never move)
    trialbeads = np.zeros(np.sum(link[:, 0]), dtype=np.uint)
    j = 0
    for i_link in xrange(nbeads-1):
        if link[i_link, 0]:
            trialbeads[j] = i_link+1
            j+=1

    print 'beadgroups = ', beadgroups[trialbeads]
    # print 'flexlink = ', flexlink
    # print 'trialbeads = ', trialbeads
    # print "cg_coor = \n", cg_coor
    # print cg_dna._coor[frame, :, 0]
    cg_dna.setCoor(cg_coor)  #; print cg_dna._coor[frame, :, 0]

    cg_dna.write_pdb("cgDNA.pdb", frame, 'w')

    vecXYZ = np.zeros((3, nbeads, 3))
    vecXYZ[0] = [1, 0, 0]
    vecXYZ[1] = [0, 1, 0]
    vecXYZ[2] = [0, 0, 1]

    #~~~~~~~~~~ END OF DNA SECTION, BEGINNING OF PROTEIN SECTION ~~~~~~~~~~~~#

    # load the protein into its own sasmol object
    aa_pro = sasmol.SasMol(0)
    error, mask = aa_all.get_subset_mask('((moltype[i] == "protein"))')
    error = aa_all.copy_molecule_using_mask(aa_pro, mask, frame)
    aa_pro.write_pdb("aaPro.pdb", frame, 'w')
    print 'protein atoms = ', aa_pro.natoms()

    # coarse-grain the proteins
    print "coarse-graining the protein, this may take awhile..."
    cg_pro = sasmol.SasMol(0)
    error, mask = aa_pro.get_subset_mask("(name[i] == 'CA')")
    error = aa_pro.copy_molecule_using_mask(cg_pro, mask, frame)
    cg_pro.write_pdb("cgPro.pdb", frame, 'w')

    # get the masks for the move groups and each group of all atom proteins
    cg_pgroup_masks = []  # masks to separate the cg_protein into NCP groups
    aa_pgroup_masks = []  # masks to separate the aa_protein into NCP groups
    aa2cg_pgroup_masks = []  # masks to get the CA atoms from the aa
    all_proteins = []

    # after using psfgen, identicle proteins end up with the same chain
    # but different segname (this creates problems for non-unique segname)
    # (pre, post, btwn) = ("((chain[i]=='", "'))", " or ")
    (pre, post, btwn) = ("((segname[i]=='", "'))", " or ")

    for group in pro_groups:
        basis_filter = ""
        for chain in group:
            basis_filter += pre + chain + post + btwn
        basis_filter = basis_filter[0:-4]  # remove the last btwn
        error, mask = cg_pro.get_subset_mask(basis_filter)
        cg_pgroup_masks.append(mask)

        aa_pro_group = sasmol.SasMol(0)
        error, mask = aa_pro.get_subset_mask(basis_filter)
        aa_pgroup_masks.append(mask)
        error = aa_pro.copy_molecule_using_mask(aa_pro_group, mask, frame)
        all_proteins.append(aa_pro_group)
        error, mask = aa_pro_group.get_subset_mask("(name[i] == 'CA')")
        aa2cg_pgroup_masks.append(mask)

    # combined the protein groups into move groups
    move_masks = []
    for i in xrange(len(cg_pgroup_masks)):
        move_masks.append(np.copy(cg_pgroup_masks[i]))
    for i in xrange(len(move_masks)):
        for j in xrange(i+1, len(move_masks)):
            move_masks[i] += move_masks[j]

    # group_indices = []
    # for group_mask in group_masks:
    #     ind = np.transpose(np.nonzero(group_mask))
    #     group_indices.append(ind[:, 0])

    # print "M (end of 'make_cgDNA_model')\n", M
    return (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, all_beads, all_proteins,
            trialbeads, beadgroups, dna_masks, aa_pgroup_masks, move_masks,
            cg_pgroup_masks)

def recover_aaPro_model(aa_pgroup_masks, cg_pro_final, cg_pro_orig, cg_pgroup_masks, all_proteins, aa_pro):
    """
    The purpose of this is to recover the all atom representation of the proteins in the new locations after rotations
    """
    test = sasmol.SasMol(0)
    error, test_mask = aa_pro.get_subset_mask("(name[i] == 'CA')")
    error = aa_pro.copy_molecule_using_mask(test, test_mask, 0)
    diffBefore = test.coor()[0] - cg_pro_orig.coor()[0]

    print 'diff before:\n', diffBefore

    for (i, mask) in enumerate(aa_pgroup_masks):
        cg_pro_group_final = sasmol.SasMol(0)                                             # define a sasmol object
        cg_pro_group_orig = sasmol.SasMol(0)
        error = cg_pro_final.copy_molecule_using_mask(cg_pro_group_final, cg_pgroup_masks[i], 0) # get the cg protein group
        error = cg_pro_orig.copy_molecule_using_mask(cg_pro_group_orig, cg_pgroup_masks[i], 0) # get the cg protein group

        coor_sub_1 = cg_pro_group_orig.coor()[0]
        com_sub_1 = cg_pro_group_orig.calccom(0)

        coor_sub_2 = cg_pro_group_final.coor()[0]
        com_sub_2 = cg_pro_group_final.calccom(0)

        all_proteins[i].align(0, coor_sub_2, com_sub_2, coor_sub_1, com_sub_1)

        error = aa_pro.set_coor_using_mask(all_proteins[i], 0, aa_pgroup_masks[i])

    error = aa_pro.copy_molecule_using_mask(test, test_mask, 0)
    diffAfter = test.coor()[0] - cg_pro_final.coor()[0]
    print 'diff after:\n', diffAfter
    print 'pause'
    return aa_pro



def recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, allBeads, masks):
    error =[]
    cg_natoms = np.copy(cg_dna.natoms())
    coor = np.copy(cg_dna.coor()[0]) #; print 'coor = ', coor
    #av = checkMag(vecXYZ) #; print 'avMag = ', av

    # split vecXYZ into three matrices
    vecX = np.copy(vecXYZ[0])
    vecY = np.copy(vecXYZ[1])
    vecZ = np.copy(vecXYZ[2])

    for i in xrange(cg_natoms):
        #s print 'allBeads[i].com (before) = \n', allBeads[i].com()

        #R is the matrix that would align the rotated coordinates back to the original
        R = align2xyz(vecX[i, :], vecY[i, :], vecZ[i, :]) #; print 'R = ', R

        #M will rotate cooridates from the original reference to the rotated coordinates of the bead then translate the com to the beads com
        M = R.transpose()
        M[3, :3] = coor[i, :] #; print 'M = ', M  # put the translation to the new com into the matrix

        r, c = allBeads[i].coor()[0].shape
        beadCoor = np.ones((r, 4))           # initialize the beadCoor matrix
        origCoor = np.ones((1, r, 3))
        origCoor = allBeads[i].coor()
        beadCoor[:, :3] = np.copy(origCoor[0]) #; print 'beadCoor (before) = \n', beadCoor
        beadCoor = np.dot(beadCoor, M) #; print 'beadCoor (after) = \n', beadCoor

        tmpCoor = np.zeros((1, r, 3))        # initialize the tmpCoor
        tmpCoor[0] = beadCoor[:, :3]      #; print 'tmpCoor = ', tmpCoor#[:, 0:3]
        allBeads[i].setCoor(tmpCoor)

        #recombine all the beads back into one pdb
        e = aa_dna.set_coor_using_mask(allBeads[i], 0, masks[i])
        error.append(e)

        #reset the bead coordinates for the next iteration
        allBeads[i].setCoor(origCoor)

    return error, aa_dna

def align2xyz(vecX, vecY, vecZ):

    tmp_coor = np.zeros((2, 4))
    tmp_coor[:, 3] = 1
    tmp_coor[1, 0:3] = vecZ
    A1 = align2z(tmp_coor)

    newX = np.dot(vecX, A1[0:3, 0:3]) #; print 'newX = ', newX
    #newY = np.dot(vecY, A1[0:3, 0:3]) ;# print 'newY = ', newY
    assert newX[2] < 10e-5, "ERROR!!! z-component of newX is not zero and it should be"

    thetaZ_x = -np.arctan2(newX[1], newX[0])  #;  print 'thetaZ_x = ', thetaZ_x
    #thetaZ_y = -np.arctan2(newY[1], -newY[0])  ;#  print 'thetaZ_y = ', thetaZ_y

    A2 = rotate4x4('z', thetaZ_x)  #; print 'A2 = ', A2

    A = np.dot(A1, A2)

    newY = np.dot(vecY, A[0:3, 0:3]) #; print 'finalY = ', newY
    assert newY[0]+newY[2] < 1+10e-5, "ERROR!!! newY is not aligned to the y-axis and it should be"

    return A

def rotate4x4(axis, theta):
    R = np.eye(4)
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

    # print 'R = ', R

    return R


def move2origin(coor4):
    T    = np.eye(4, dtype=np.float)
    Ti   = np.copy(T)
    #s print 'coor passed to T, Ti:\n', coor4
    #s print 'size = ', coor4.shape
    if not coor4.shape < (1, 4):
        T[3, 0:3] = -coor4[0, 0:3]
        Ti[3, 0:3] = coor4[0, 0:3]

    #s print 'T:\n', T
    #s print 'test:\n', np.dot(coor4, T)
    #s print 'Ti:\n', Ti
    #s print 'testI:\n', np.dot(np.dot(coor4, T), Ti)
    return (T, Ti)

def align2z(coor4):
    '''
    function designed to align the axis connecting the first 2 coordinates of an array
    of (1, 4) coodinate vectors to the z-axis
    '''
    A = np.eye(4, dtype=np.float)
    #s Axz  = np.eye(4, dtype=np.float)
    #s Az   = np.eye(4, dtype=np.float)
    assert all(coor4[0] == [0., 0., 0., 1., ]), "coordinates passed to align2z were not translated to the origin"

    if coor4.shape > (1, 4):
        (u, v, w) = coor4[1, 0:3]
        small = 1E-14 # small limit to make sure not dividing by zero

        d1 = np.sqrt(u**2+v**2)
        d2 = np.sqrt(u**2+v**2+w**2)
        if d1 > small:
            (A[0][0], A[0][1], A[0][2]) = (u/d1*w/d2, -v/d1, u/d2)
            (A[1][0], A[1][1], A[1][2]) = (v/d1*w/d2, u/d1, v/d2)
            (A[2][0], A[2][1], A[2][2]) = (-d1/d2, 0, w/d2)   #; print 'A1:\n', A

        #s # align to the x-z plane
        #s if v > small:           # check if already aligned to xz-plane
        #s         (Axz[0][0], Axz[0][1]) = (u/d1, -v/d1)
        #s         (Axz[1][0], Axz[1][1]) = (v/d1, u/d1)  #; print 'Axz= \n', Axz
        #s
        #s # align to the z-axis
        #s if d1 > small:
        #s         (Az[0][0], Az[0][2]) = (w/d2, d1/d2)
        #s         (Az[2][0], Az[2][2]) = (-d1/d2, w/d2)  #; print 'Az= \n', Az
    else:
        print 'no point to align'
        #s A = np.dot(Axz, Az) #;   print 'A3=\n', A
    return A


def beadRotate(coor3, vecXYZ, thetas, nSoft, p_coor3):
    '''
    this function is designed to generate a modified version of the input coordinates (coor3)
    it moves the all coordinates so the first is at the origin
    aligns the first two coordinates to be along the z-axis
    performs the rotation
    it then undoes the alignment to the z-axis and translates all coordinates returning the first coordinate
    to where it started
    '''
    # debug input/output
    #print 'vecXYZ = (before mod, before assign)\n', vecXYZ

    # add a fourth column of ones to each bead origin and orientation vectors (this incorporates translations)
    (natoms, col) = coor3.shape

        # initialize vector arrays for coordinates and orientation vectors making  them into 4 component vectors
    coor4 = np.ones((natoms, 4), np.float)
    X = np.copy(coor4)
    Y = np.copy(coor4)
    Z = np.copy(coor4)
    coor4[:, 0:3] = coor3     #; print 'coor4 = ', coor4
    X[:, 0:3] = np.copy(vecXYZ[0])  #;print 'X = \n', X
    Y[:, 0:3] = np.copy(vecXYZ[1])
    Z[:, 0:3] = np.copy(vecXYZ[2])

    p_coor4 = np.ones((len(p_coor3), 4), np.float)
    p_coor4[:, 0:3] = p_coor3

    # create the translation-rotation matrix
    # This is intended to be multiplied from the right (unlike standard matrix multiplication) so as not to require transposing the coordinate vectors.
    # print '(tx, ty, tz) in degrees:', thetas
    tx = thetas[0]*(np.pi/180.0)
    ty = thetas[1]*(np.pi/180.0)
    tz = thetas[2]*(np.pi/180.0)
    # print '(tx, ty, tz) in radians:', (tx, ty, tz)
    cx = np.cos(tx)
    sx = np.sin(tx)
    cy = np.cos(ty)
    sy = np.sin(ty)
    cz = np.cos(tz)
    sz = np.sin(tz)

    # initialize the rotation
    # consolidated method of defining the rotation matrices
    R = np.eye(4, dtype=np.float)
    (R[0][0], R[0][1], R[0][2]) = ( cy*cz, cy*sz, -sy )
    (R[1][0], R[1][1], R[1][2]) = ( sx*sy*cz-cx*sz, sx*sy*sz+cx*cz, sx*cy )
    (R[2][0], R[2][1], R[2][2]) = ( sx*sz+cx*sy*cz, cx*sy*sz-sx*cz, cx*cy )  #; print 'R:\n', R

    # verbose method of defining the rotation matrices (slower)
    #s Rx = np.eye(4, dtype=np.float)
    #s Ry = np.eye(4, dtype=np.float)
    #s Rz = np.eye(4, dtype=np.float)
    #s
    #s (Rx[1][1], Rx[1][2]) = ( cx, sx)
    #s (Rx[2][1], Rx[2][2]) = (-sx, cx)  #; print 'Rx:\n', Rx
    #s
    #s (Ry[0][0], Ry[0][2]) = (cy, -sy)
    #s (Ry[2][0], Ry[2][2]) = (sy, cy)  #; print 'Ry:\n', Ry
    #s
    #s (Rz[0][0], Rz[0][1]) = ( cz, sz)
    #s (Rz[1][0], Rz[1][1]) = (-sz, cz)  #; print 'Rz:\n', Rz
    #s
    #s R = np.dot(np.dot(Rx, Ry), Rz) #; print "Rxyz = \n", Rxyz

    # print 'original coor:\n', coor4
    (T0, Ti0) = move2origin(coor4)  # Ti0 is NOT the transpose of T0, rather the off diag elements are negative
    coor4 = np.dot(coor4, T0)  #; print 'moved to origin coor:\n', coor4
    A = align2z(coor4)
    AR = np.dot(A, R)

    # first of nSoft rotations
    coor4 = np.dot(coor4, AR)         #; print 'aligned coor to z-axis:\n', coor4 #; print 'step 0 rotated about first angle coor:\n', coor4

    # the coarse grained beads local coordinates should not be translated, only rotated
    X = np.dot(X, AR)
    Y = np.dot(Y, AR)
    Z = np.dot(Z, AR)

    p_coor4 = np.dot(np.dot(p_coor4, T0), AR)

    # remaining of nSoft rotations
    for i in xrange(1, nSoft):
        (T, Ti) = move2origin(coor4[i:])
        coor4[i:] = np.dot(coor4[i:], T)  # move to origin
        # print 'moved2origin coor:\n', coor4
        p_coor4 = np.dot(p_coor4, T)

        coor4[i:] = np.dot(coor4[i:], R)
        # print 'step %d' %i, 'rotated coor:\n', coor4
        p_coor4 = np.dot(p_coor4, R)

        coor4[i:] = np.dot(coor4[i:], Ti) # return to original position
        # print 'returned from origin coor:\n', coor4
        p_coor4 = np.dot(p_coor4, Ti)

        # the coarse grained beads local coordinates should not be translated, only rotated
        X[i:] = np.dot(X[i:], R)
        Y[i:] = np.dot(Y[i:], R)
        Z[i:] = np.dot(Z[i:], R)


    coor4 = np.dot(coor4, A.transpose())
    coor4 = np.dot(coor4, Ti0)        # print 'returned from origin coor:\n', coor4
    X = np.dot(X, A.transpose())
    Y = np.dot(Y, A.transpose())
    Z = np.dot(Z, A.transpose())

    p_coor4 = np.dot(p_coor4, A.transpose())
    p_coor4 = np.dot(p_coor4, Ti0)

    vecXYZ[0, 1:] = X[1:, 0:3]
    vecXYZ[1, 1:] = Y[1:, 0:3]
    vecXYZ[2, 1:] = Z[1:, 0:3]



    return (coor4[1:, 0:3], vecXYZ[:, 1:], p_coor4[:, 0:3])           # this returns the modified positions and orientations for all but the first (reference) bead

def checkU(coor):
    '''
    module to generate and check the u-vectors pointing from bead to bead
    makes sure they are all the same length
    returns the u-vectors and their magnitude (l)
    '''
    erLim = 1e-3
    (r, c) = coor.shape
    u = coor[1:, :]-coor[:r-1, :]                 # u-vectors pointing from bead to bead
    lu = np.sqrt(u[:, 0]**2+u[:, 1]**2+u[:, 2]**2) # magnitues of u-vectors -> distance btwn beads
    l = np.mean(lu)                          # average distance btwn bead

    #print '\n\ncoor:\n', coor
    #print 'u-vectors:\n', u
    #print 'magnitude of u-vectors:\n', lu
    #print 'average magnitude of u =', l
    #print '\n\n'

    test = np.abs(lu-l) > erLim  # check if any of the l-lengths are different
    #s if test.any():
    #s         print 'ERROR: the beads are not uniformly spaced'
    #s         print 'u = \n', u
    #s         print test
    #s         # print 'difference = ', u[:, 2] - l
    #s         print 'distance = ', lu

    return (u, l)

def checkMag(vec):
    '''
    module to check that the magnitude of all vectors is the same
    '''
    erLim = 1e-3
    (r, c) = vec.shape
    m = np.zeros(r)
    for i in xrange(r):
        m[i] = mag(vec[i])

    avMag = np.mean(m)
    test = np.abs(m-avMag) > erLim  # check if any of the l-lengths are different
    if test.any():
        print 'ERROR: the vectors do not have uniformly magnitude'
        print 'm = \n', m
        print test
        # print 'difference = ', u[:, 2] - l
        print 'vec = ', vec

    return (avMag)

def energyBend(lpl, u, l):
    '''
    this function finds the E/kT for N beads connected by N-1 rods
    lpl is the persistence length divided by the rod length: lp/l
    u contains the N-1 unit vectors representing each rod
    '''

    (nu, col) = u.shape          # nu = nbeads - 1
    uk0 = u[:nu-1, :]/l          #; print 'u_k= \n', uk0
    uk1 = u[1:, :]/l             #; print 'u_k= \n', uk1
    #s print 'uk0*uk1\n:', uk0*uk1

    # checkMag(u)       # for debugging

    res = (nu-1) - np.sum(uk0*uk1)  # calculate the sum of their products

    return res*lpl

def energyWCA(w, coor, wca0, trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    '''
    wca1 = np.copy(wca0)
    test = 2.**(1./6.)*w
    (N, col) = coor.shape
    for i in xrange(trial_bead, N):
        ri = coor[i, :]
        for j in xrange(0, i):
            rj = coor[j, :]
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)
            #s print '(2^(1/6)*w, rij) = ', (test, rij)
            if rij < test:
                wca1[i, j] = (w/rij)**12.-(w/rij)**6.+0.25
    res = 4*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def energyWCA_mixed(w_d, coor_d, wca0_d, trial_bead, w_p, coor_p, wca0_p, group_mask_p, wca0_dp):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    '''
    w_dp = np.mean([w_d, w_p])

    wca1_d = np.copy(wca0_d)
    wca1_dp = np.copy(wca0_dp)
    wca1_p = np.copy(wca0_p)

    test_d = 2.**(1./6.)*w_d
    test_p = 2.**(1./6.)*w_p
    test_dp = 2.**(1./6.)*w_dp

    (N_d, col) = coor_d.shape
    (N_p, col) = coor_p.shape

    # calculate the distance between moved DNA beads and all other DNA beads
    for i in xrange(trial_bead, N_d):    # only calculate the distances that changed
        ri = coor_d[i, :]                # get the coordinates for the ith bead
        for j in xrange(i):
            rj = coor_d[j, :]            # get the coordinates for the jth bead
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith and jth beads
            if rij < test_d:
                wca1_d[i, j] = (w_d/rij)**12.-(w_d/rij)**6.+0.25

    # calculate the distance between DNA beads and proteins
    for i in xrange(N_d):
        ri = coor_d[i, :]
        for j in xrange(N_p):
            rj = coor_p[j, :]
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
            if rij < test_dp:
                wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25

    # calculate the distance between proteins and proteins
    for i in xrange(N_p):
        ri = coor_p[i, :]
        for j in xrange(i):
            rj = coor_p[j, :]
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
            if rij < test_p:
                wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25

    res_d = 4 * np.sum(wca1_d)
    res_p = 4 * np.sum(wca1_p)
    res_dp = 4 * np.sum(wca1_dp)

    res = res_d + res_p + res_dp
    return (res, wca1_d, wca1_p, wca1_dp)

def energyWCA_mixedFaster(w_d, coor_d, wca0_d, trial_bead, w_p, coor_p, wca0_p, group_mask_p, wca0_dp, initial=False):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads

    this method skips unnecessary distances, compared to other method and difference was attributable to float error
    '''
    w_dp = np.mean([w_d, w_p])

    wca1_d = np.copy(wca0_d)
    wca1_dp = np.copy(wca0_dp)
    wca1_p = np.copy(wca0_p)

    test_d = 2.**(1./6.)*w_d
    test_p = 2.**(1./6.)*w_p
    test_dp = 2.**(1./6.)*w_dp

    (N_d, col) = coor_d.shape
    (N_p, col) = coor_p.shape

    # calculate the distance between moved DNA beads and all other DNA beads
    for i in xrange(trial_bead, N_d):    # only calculate the distances that changed
        ri = coor_d[i, :]                # get the coordinates for the ith bead
        for j in xrange(i):
            rj = coor_d[j, :]            # get the coordinates for the jth bead
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith and jth beads
            if rij < test_d:
                wca1_d[i, j] = (w_d/rij)**12.-(w_d/rij)**6.+0.25

    # get the indices of the moved and stationary proteins
    ind_moved_p = mask2ind(group_mask_p)
    ind_stationary_p = mask2ind(-(group_mask_p-1))

    if initial:
        for i in xrange(N_d):           # want all the DNA
            ri = coor_d[i, :]           # get the coordinates for the ith bead
            for j in xrange(N_p):       # want all the protein
                rj = coor_p[j, :]
                rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
                if rij < test_dp:
                    wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25
        for i in xrange(N_p):
            ri = coor_p[i, :]
            for j in xrange(i):
                rj = coor_p[j, :]
                rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
                if rij < test_p:
                    wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25
    else:
        # calculate the distance between moved DNA and stationary proteins
        if len(ind_stationary_p) > 0:
            for i in xrange(trial_bead, N_d):    # only want the DNA that moved
                ri = coor_d[i, :]                # get the coordinates for the ith bead
                for j in ind_stationary_p:       # only want the protein that was stationary
                    rj = coor_p[j, :]
                    rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
                    if rij < test_dp:
                        wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25
        # calculate the distance between stationary DNA and moved proteins
        if len(ind_moved_p) > 0:
            for i in xrange(trial_bead):         # only want the DNA that did not move
                ri = coor_d[i, :]                # get the coordinates for the ith bead
                for j in ind_moved_p:            # only want the protein that was stationary
                    rj = coor_p[j, :]
                    rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
                    if rij < test_dp:
                        wca1_dp[i, j] = (w_dp/rij)**12.-(w_dp/rij)**6.+0.25

        # calculate the distance between moved proteins and stationary proteins
        if len(ind_moved_p) > 0 and len(ind_stationary_p) > 0:
            for i in ind_moved_p:
                ri = coor_p[i, :]
                for j in ind_stationary_p:
                    rj = coor_p[j, :]
                    rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)  # calculate the distance between the ith DNA and jth protein beads
                    if rij < test_p:
                        wca1_p[i, j] = (w_p/rij)**12.-(w_p/rij)**6.+0.25

    res_p = 4 * np.sum(wca1_p)
    res_d = 4 * np.sum(wca1_d)
    res_dp = 4 * np.sum(wca1_dp)

    res = res_d + res_p + res_dp
    return (res, wca1_d, wca1_p, wca1_dp)

def FenergyWCA_mixed(w_d, coor_d, wca0_d, trial_bead0, w_p, coor_p, wca0_p, group_mask_p, wca0_dp, initial=False):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN

    this method skips unnecessary distances, compared to other method and difference was attributable to float error
    '''
    #import sys ; sys.path.append('./')  #import sys ; sys.path.append('/home/schowell/Dropbox/gw_phd/code/pylib/sassie/')
    #import electrostatics

    w_dp = np.mean([w_d, w_p])
    wca1_d = np.copy(wca0_d)
    wca1_dp = np.copy(wca0_dp)
    wca1_p = np.copy(wca0_p)

    (N_d, col) = coor_d.shape
    (N_p, col) = coor_p.shape

    # get the indices of the moved and stationary proteins
    ind_moved_p = mask2ind(group_mask_p)
    ind_stationary_p = mask2ind(-(group_mask_p-1))
    if len(ind_moved_p):
        imoved0 = np.min(ind_moved_p)
    else:
        imoved0 = len(ind_stationary_p) - 1
    imoved1 = imoved0 + 1
    trial_bead1 = trial_bead0 + 1

    # calculate the distance between moved DNA beads and all other DNA beads
    wca1_d = collision.wca_d(coor_d, trial_bead1, w_d, wca1_d)    #  increment trial bead by one because it will index from 1 instead of 0
    res_d = 4.*np.sum(wca1_d)

    if initial:
        wca1_p = collision.wca_d(coor_p, 1, w_p, wca1_p)
        res_p = 4 * np.sum(wca1_p)

        wca1_dp = collision.wca_nbym(coor_d, coor_p, 1, N_d, 1, N_p, w_dp, wca1_dp)
        res_dp = 4 * np.sum(wca1_dp)

    else:
        # calculate the distance between moved proteins and stationary proteins
        if len(ind_moved_p) > 0 and len(ind_stationary_p) > 0:
            wca1_p = collision.wca_nbyn(coor_p, imoved1, N_p, 1, imoved0, w_p, wca1_p)    #  increment trial bead by one because it will index from 1 instead of 0
        res_p = 4 * np.sum(wca1_p)

        # calculate the distance between moved DNA and stationary proteins
        if len(ind_stationary_p) > 0:
            wca1_dp = collision.wca_nbym(coor_d, coor_p, trial_bead1, N_d, 1, imoved0, w_dp, wca1_dp)
        # calculate the distance between stationary DNA and moved proteins
        if len(ind_moved_p) > 0:
            wca1_dp = collision.wca_nbym(coor_d, coor_p, 0, trial_bead0, imoved1, N_p, w_dp, wca1_dp)
        res_dp = 4 * np.sum(wca1_dp)

    #s res = res_d + res_p + res_dp # I think the dna-protein calculation is throwing things off
    res = res_d + res_p
    return (res, wca1_d, wca1_p, wca1_dp)

def FenergyWCA(w, coor, wca0, trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN
    '''
    import sys ; sys.path.append('./')  #import sys ; sys.path.append('/home/schowell/Dropbox/gw_phd/code/pylib/sassie/')
    import electrostatics

    wca1 = np.copy(wca0)
    wca1 = electrostatics.calcwca(coor, trial_bead+1, w, wca1)
    #  increment trial bead by one because it will index from 1 instead of 0

    res = 4.*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def mask2ind(mask):
    '''
    convert a mask to an array of indicies
    '''

    mask_indices = np.nonzero(mask*np.arange(1, len(mask)+1))[0]

    return mask_indices

def test_wca(val1, val2, out=False):
    if val1.shape != val2.shape:
        print 'ERROR: val1 and val2 are not comparable'
        return

    diff = (val1 - val2)
    per_diff = diff / val1
    if np.any(np.isnan(per_diff)):
        per_diff[np.isnan(per_diff)] = 0
    sum_diff = np.sum(per_diff)
    max_diff = np.max(np.abs(per_diff))

    print 'val1:\n', val1
    print 'val2:\n', val2
    print 'difference:\n', diff
    print 'percent difference:\n', per_diff
    print 'total percent difference:\n', sum_diff
    print 'maximum percent difference:\n', max_diff

    if out:
        return (per_diff, max_diff, diff)

    return ''

def dna_mc(nsteps, cg_dna, cg_pro, vecXYZ, lp, w, theta_max, trialbeads, beadgroups, group_masks, nSoft=3, f=True):
    #def dna_mc(nsteps, cg_dna, vecXYZ, lp, w, theta_max, trialbeads, nSoft=3, f=True):
    '''
    this function perform nsteps Monte-Carlo moves on the cg_dna
    '''
    dcdOutFile_d = cg_dna.open_dcd_write("cg_dna_moves.dcd")
    dcdOutFile_p = cg_pro.open_dcd_write("cg_pro_moves.dcd")

    nbeads = cg_dna.natoms()
    nbeads_p = cg_pro.natoms()
    nflex = trialbeads.size
    xyz = np.copy(vecXYZ)
    coor = np.copy(cg_dna.coor()[0])
    p_coor = np.copy(cg_pro.coor()[0])

    (u, l) = checkU(coor) # get the vectors pointing from bead to bead, u, and the average distance, l
    lpl = lp/l  # setup the presistence length paramater
    w_p = w / 20.  # this is arbitrary, need to look learn more about proteins

    # calculate the energy of the starting positions
    wca0 = np.zeros((nbeads, nbeads))
    wca0_p = np.zeros((nbeads_p, nbeads_p))
    wca0_dp = np.zeros((nbeads, nbeads_p))

    Ub0 = energyBend(lpl, u, l)
    if f:
        # (Uwca0, wca0) = FenergyWCA_mixed(w, coor, wca0, 0, cg_pro.coor()[0])
        (Uwca0, wca0, wca0_p, wca0_dp) = FenergyWCA_mixed(w, coor, wca0, 0, w_p, p_coor, wca0_p, group_masks[1], wca0_dp, True)
    else:
        (Uwca0, wca0, wca0_p, wca0_dp) = energyWCA_mixed(w, coor, wca0, 0, w_p, p_coor, wca0_p, group_masks[1], wca0_dp)

    U_T0 = Ub0 + Uwca0
    #s wca1 = np.copy(wca0)

    a = 0 # times configuration was accepted
    r = 0 # times configuration was rejected

    for i in xrange(nsteps):

        trial_bead = trialbeads[int((nflex)*random.random())]  #; print 'trial_bead =', trial_bead
        # trial_bead = trialbeads[i%len(trialbeads)]

        thetaZ_max = np.float(theta_max) / 50. # this should be something small
        thetaX = theta_max * random.random() - theta_max/2
        thetaY = theta_max * random.random() - theta_max/2
        thetaZ = thetaZ_max * random.random() - thetaZ_max/2
        thetaXYZ = [thetaX/nSoft, thetaY/nSoft, thetaZ/nSoft]  #; print 'thetaXYZ = ', thetaXYZ

        # get the protein coordinates for the group that will be transformed
        if beadgroups[trial_bead] == len(group_masks):
            p_coor_rot = np.zeros((1, 1, 3))
            p_mask = np.zeros(len(group_masks[0]))
        else:
            p_mask = group_masks[beadgroups[trial_bead]]
            p_indices = mask2ind(p_mask)
            p_coor_rot = p_coor[p_indices]

        (coor[trial_bead:], xyz[:, trial_bead:], p_coor_rot) = beadRotate(coor[trial_bead-1:], xyz[:, trial_bead-1:], thetaXYZ, nSoft, p_coor_rot) # generate a newly rotated model

        # store the rotated protein coordinates
        if beadgroups[trial_bead] < len(group_masks):
            p_coor[p_indices] = p_coor_rot

        # calculate the change in energy (dU) and boltzman factor (p) for the new model
        (u, l) = checkU(coor)
        Ub1 = energyBend(lpl, u, l)
        if f:
            (Uwca1, wca1, wca1_p, wca1_dp) = FenergyWCA_mixed(w, coor, wca0, trial_bead, w_p, p_coor, wca0_p, p_mask, wca0_dp, True)
            # (Uwca1, wca1) = FenergyWCA(w, coor, wca0, trial_bead) # just DNA
        else:
            (Uwca1, wca1, wca1_p, wca1_dp) = energyWCA_mixedFaster(w, coor, wca0, trial_bead, w_p, p_coor, wca0_p, p_mask, wca0_dp)  # DNA-protein
            # (Uwca1, wca1) = energyWCA(w, coor, wca0, trial_bead) # just DNA
            # tic = time.time()
            # (Uwca1, wca1, wca1_p, wca1_dp) = energyWCA_mixed(w, coor, wca0, trial_bead, w_p, p_coor, wca0_p, p_mask, wca0_dp)
            # toc = time.time() - tic ; print 'WCA long time =', toc, 'seconds'

        U_T1 =  Ub1 + Uwca1
        dU = U_T1 - U_T0
        p = np.exp(-dU)

        test = random.random()

        # if accepted write new coordinates, else write old again
        # if True:
        if test < p:
            # print 'wrote new dcd frame (end of loop', i, ' trial_bead=', trial_bead, ')'+' accepted new configuration\n'
            a += 1                                     # increment the accepted counter
            cg_dna.coor()[0] = np.copy(coor)           # store the dna coordinates
            cg_pro.coor()[0] = np.copy(p_coor)
            vecXYZ = np.copy(xyz)                      # store the dna orientations
            wca0 = np.copy(wca1)                       # update the WCA energy for the DNA
            wca0_dp = np.copy(wca1_dp)                 # update the WCA energy for the DNA-protein
            wca0_p = np.copy(wca1_p)                   # update the WCA energy for the protein
            U_T0 = U_T1                                # update the total energy
            cg_dna.write_dcd_step(dcdOutFile_d, 0, 0)
            cg_pro.write_dcd_step(dcdOutFile_p, 0, 0)
        else :
            # print 'wrote new dcd frame (end of loop', i, ' trial_bead=', trial_bead, ')'+' rejected new configuration\n'
            r += 1                                     # increment the rejected counter
            coor = np.copy(cg_dna.coor()[0])           # reset the dna coordinates
            p_coor = np.copy(cg_pro.coor()[0])         # reset the protein coordinates
            xyz = np.copy(vecXYZ)                      # reset the dna orientations

        # cg_dna.write_dcd_step(dcdOutFile_d, 0, 0)
        # cg_pro.write_dcd_step(dcdOutFile_p, 0, 0)
        # print 'coor = \n', coor

    # close output files
    cg_dna.close_dcd_write(dcdOutFile_d)
    cg_pro.close_dcd_write(dcdOutFile_p)
    #s outData.close()

    return (cg_dna, vecXYZ, cg_pro, a, r)

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

    longCoor = np.zeros((1, natoms, 3), np.float)    # initialize the long DNA coordinates
    longCoor[0][:, 2] = range(natoms)      # set the z-values to the index of the array
    longCoor *= l                   # scale the z-values to the right seperation
    # print longCoor[-5:]

    longDNA.setCoor(longCoor)
    longDNA.setElement(['C']*natoms)

    vecXYZ = np.zeros((3, natoms, 3))
    vecXYZ[0] = [1, 0, 0]
    vecXYZ[1] = [0, 1, 0]
    vecXYZ[2] = [0, 0, 1]
    # n = L/l                         # number of times need to repeat the grain
    # print '(l, L, n)', (l, L, n)

    return (longDNA, vecXYZ)

def mag(vec):

    (r, ) = vec.shape
    sumSquares = 0.0
    for i in xrange(r):
        sumSquares += vec[i]**2

    return np.sqrt(sumSquares)

def closeAll():

    for x in xrange(100):
        plt.close()

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-i", "--iters", default="100", type=int, help="number of times to repeat n-steps (program saves coordinates after each iteration)")
    parser.add_argument("-n", "--nsteps", default="1", type=int, help="number of Monte Carlo steps in each iteration (program saves coordinates after this many steps)")
    parser.add_argument("-s", "--show", action="store_true", help="show a plot of the configuration after each iteration")
    parser.add_argument("-tm", "--theta_max", default="10", type=np.float, help="max rotation angle for the x and y rotation, max z rotation is theta_max/50")
    parser.add_argument("-ns", "--nsoft", default="2", type=int, help="number of bead over which to apply any rotation (i.e. soften the bending)")
    parser.add_argument("-bp", "--bp_perBead", default="5", type=int, help="number of bp for each coarse-grained bead")
    parser.add_argument("-c", "--cgpdb", help="coarse-grained pdb file")

    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument("-f", "--for_collide", action="store_true", help="use fortran module to calculate the atom overlap")
    group1.add_argument("-py", "--py_collide", action="store_true", help="use python module to calculate the atom overlap")

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-l", "--Llp", type=int, help="length of desired coarse-grain long DNA to simulate (in units of persistence lengths)")
    group2.add_argument("-p", "--pdb", default="dna.pdb", help="all atom pdb file")

    # parser.add_argument("--option", metavar='', help='')
    return parser.parse_args()



def main():

    iters = ARGS.iters
    nsteps = ARGS.nsteps
    theta_max = ARGS.theta_max
    bp_perBead = ARGS.bp_perBead
    nSoft = ARGS.nsoft
    lp = 530      # persistence length  (lp = 530A)
    w = 46        # width of chain in A (w = 46A)l

    if ARGS.py_collide:
        f = False
        print 'using python collision calculation'
    else:
        f = True
        print 'using fortran collision calculation'


    if ARGS.Llp:
        print 'simulating long DNA'
        Llp = ARGS.Llp    # Llp is the length in units of lp: Llp = L/lp
        L = Llp*lp  #
        (cg_dna, vecXYZ) = makeLongDNA(Llp) # use this to make long cgDNA
        all_atom_pdb = '%d_llp' % Llp
    elif ARGS.pdb:
        print 'loading pdb'
        all_atom_pdb = ARGS.pdb

        if ARGS.pdb == '1zbb_tetra_half.pdb':
            dna_chains = ['I', 'J']
            dna_resids.append([1, 347])
            dna_resids.append([347, 1])
            # continuous flexible residues, recall excludes upper lim: [a, b)
            l = [range(1, 31), range(165, 196), range(333, 348)]
            # proteing groups
            pro_groups = []
            pro_groups.append(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
            pro_groups.append(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])

        elif  ARGS.pdb ==  'c11_withTails.pdb':
            '''The C11 nucleosome tetramer with all the protein tails'''
            dna_resids = ['J', 'L']
            dna_resids.append([1, 693])
            dna_resids.append([693, 1])
            # continuous flexible residues of DNA on first chain
            # recall that range(a, b) excludes upper lim: [a, b)
            l = [range(1, 31), range(167, 198), range(334, 365),
                 range(501, 532), range(667, 694)]
            pro_groups =  []
            pro_groups.append(['m0', 'n0', 'o0', 'p0',
                               'q0', 'r0', 's0', 't0'])
            pro_groups.append(['M1', 'N1', 'O1', 'P1',
                               'Q1', 'R1', 'S1', 'T1'])
            pro_groups.append(['A1', 'B1', 'C1', 'D1',
                               'E1', 'F1', 'G1', 'H1'])
            pro_groups.append(['a0', 'b0', 'c0', 'd0',
                               'e0', 'f0', 'g0', 'h0'])

        elif ARGS.pdb == '1zbb_tetra_corrected.pdb':
            dna_chains = ['L', 'J']
            dna_resids.append([1, 694])
            dna_resids.append([694, 1])
            # continuous flexible residues of DNA on first chain
            # recall that range(a, b) excludes upper lim: [a, b)
            l = [range(1, 31), range(167, 198), range(334, 365),
                 range(501, 532), range(667, 695)]
            # proteing groups
            pro_groups = []
            pro_groups.append(['m', 'n', 'o', 'p', 'q', 'r', 's', 't'])
            pro_groups.append(['M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'])
            pro_groups.append(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
            pro_groups.append(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])

        elif ARGS.pdb == 'dna.pdb':
            #dummy dna file
            dna_chain1 = 'A'
            dna_chain2 = 'B'
            dna_resid1 = np.array([1, 60])
            dna_resid2 = np.array([120, 61])
            l = [range(1, 120)]
            pro_groups = []
        else:
            print "\n>>> ERROR, new file, need to set chain1, chain2, resid1, resid2, and l <<<\n"

        tic = time.time()

        (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, all_beads, all_proteins,
         trialbeads, beadgroups, dbead_masks, aa_pgroup_masks, move_masks,
         cg_pgroup_masks) = make_cg_model(all_atom_pdb, dna_chains, dna_resids,
                                          l, pro_groups, bp_perBead)

        cg_pro_orig = sasmol.SasMol(0)
        error = cg_pro.copy_molecule_using_mask(cg_pro_orig, np.ones(len(cg_pro.coor()[0])), 0)
        toc = time.time() - tic ; print 'coarse-grain time =', toc, 'seconds'

    else:
        print " ERROR! don't know what to do"

    # for bead in trialbeads:
    #    print 'bead', bead, 'group', beadgroups[bead-1]

    (nbeads, c) = cg_dna.coor()[0].shape
    # print 'nbeads = ', nbeads
    coorCopy = np.ones((nbeads, 4))
    coorCopy[:, :3] = np.copy(cg_dna.coor()[0])
    #s print cg_dna.coor()[0]

    rg_lp = np.zeros(iters)
    re_lp = np.zeros(iters)
    rg0 = cg_dna.calcrg(0)/lp
    re0 = mag(cg_dna.coor()[0, -1]-cg_dna.coor()[0, 0])/lp

    if ARGS.show:
        import pylab
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        # ax = fig.add_subplot(iters, 1, 1, projection='3d')
        ax = Axes3D(fig)
        ax.scatter(cg_dna.coor()[0, :, 0], cg_dna.coor()[0, :, 1], cg_dna.coor()[0, :, 2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fig.suptitle('After %d steps Rg/lp=%f  Re/lp=%f' %(0, rg0, re0))
        plt.show()

    print 'iter:', 0, 'of', iters, ' (a, r) = (  , )   rg/lp=', rg0, 're/lp=', re0

    # setup output
    timestr = time.strftime("%y%m%d_%H%M%S")
    fName = 'output/' + timestr + '_dnaMoves.o'
    outData = open(fName, 'a')
    outData.write("# pdb=%s\t iters=%d\t nsteps=%d\t nSoft=%d\t theta_max=%f \n# moves\t rg/lp\t\t re/lp\t\t a\t r\n" %(all_atom_pdb, iters, nsteps, nSoft, theta_max))

    dna_dcd_name = 'output/' + timestr + '_aaDNA_Moves.dcd'
    aa_dna_dcdOutFile = aa_dna.open_dcd_write(dna_dcd_name)
    pro_dcd_name = 'output/' + timestr + '_aaPro_Moves.dcd'
    aa_pro_dcdOutFile = aa_pro.open_dcd_write(pro_dcd_name)

    tic = time.time()
    for i in xrange(iters):
        (cg_dna, vecXYZ, cg_pro, a, r) = dna_mc(
            nsteps, cg_dna, cg_pro, vecXYZ, lp, w, theta_max, trialbeads,
            beadgroups, move_masks, nSoft, f)
        print 'cg_pro.coor()[0]:\n',  cg_pro.coor()[0]
        print 'cg_pro_orig.coor()[0]:\n',  cg_pro_orig.coor()[0]

        rg_lp[i] = cg_dna.calcrg(0)/lp
        re_lp[i] = mag(cg_dna.coor()[0, -1]-cg_dna.coor()[0, 0])/lp
        print 'iter:', i+1, 'of', iters, ' (a, r) = ', (a, r), 'rg/lp=', rg_lp[i], 're/lp=', re_lp[i]

        outData.write("%d\t%f\t%f\t%d\t%r\n" %((i+1)*nsteps, rg_lp[i], re_lp[i], a, r) )                #output

        if ARGS.show:
            fig = pylab.figure()
            ax = Axes3D(fig)
            ax.scatter(cg_dna.coor()[0, :, 0], cg_dna.coor()[0, :, 1], cg_dna.coor()[0, :, 2])
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            fig.suptitle('After %d steps Rg/lp=%f  Re/lp=%f' %((i+1)*nsteps, rg_lp[i], re_lp[i]))

        # recover an all atom representation
        # DNA
        error, aa_dna = recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, all_beads, dbead_masks)
        aa_dna.write_dcd_step(aa_dna_dcdOutFile, 0, 0)
        # Protein
        for (i, mask) in enumerate(aa_pgroup_masks):
            # define a temporary sasmol object for the final cg protein group
            cg_pro_group_final = sasmol.SasMol(0)
            error = cg_pro.copy_molecule_using_mask(cg_pro_group_final,
                                                    cg_pgroup_masks[i], 0)

            # get the sub_1 inputs for the align function
            com_sub_1 = cg_pro_group_final.calccom(0)
            cg_pro_group_final.center(0)
            coor_sub_1 = cg_pro_group_final.coor()[0]

            all_proteins[i].center(0)
            # define a temporary sasmol object for the original cg protein group
            cg_pro_group_orig = sasmol.SasMol(0)
            error, aa2cg_mask = all_proteins[i].get_subset_mask("(name[i] == 'CA')")
            error = all_proteins[i].copy_molecule_using_mask(cg_pro_group_orig,
                                                             aa2cg_mask, 0)
            # get the sub_2 inputs for the align function
            com_sub_2 = cg_pro_group_orig.calccom(0)
            cg_pro_group_orig.center(0)
            coor_sub_2 = cg_pro_group_orig.coor()[0]

            all_proteins[i].align(0, coor_sub_2, com_sub_2,
                                  coor_sub_1, com_sub_1)
            error = aa_pro.set_coor_using_mask(all_proteins[i], 0, mask)
        #aa_pro.write_dcd_step(aa_pro_dcdOutFile, 0, 0)

        # aa_pro = recover_aaPro_model(aa_pgroup_masks, cg_pro, cg_pro_orig, cg_pgroup_masks, all_proteins, aa_pro) # <--- should use this
        aa_pro.write_dcd_step(aa_pro_dcdOutFile, 0, 0)

        os.system("cp cg_dna_moves.dcd output/" + timestr + '_cgDNA_Moves.dcd')
        os.system("cp cg_pro_moves.dcd output/" + timestr + '_cgPro_Moves.dcd')

    toc = time.time() - tic ; print 'run time =', toc, 'seconds'

    outData.close()
    # cg_dna.close_dcd_write(dcdOutFile)

    #recover an all atom representation
    # error, aa_dna = recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, allBeads, masks)
    # aa_dna.write_pdb("finalDNA.pdb", 0, 'w')

if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()
    main()