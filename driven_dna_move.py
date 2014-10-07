#!/usr/bin/env python
#!/share/apps/bin/python
# Author:   --<Steven Howell>
# Purpose:  Generate modified DNA or DNA-protein structures
# Created: 12/01/2013
# $Id: driven_dna_move.py,v 1.1 2014-10-07 20:36:48 schowell Exp $

#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sassie.sasmol.sasmol as sasmol, numpy as np
import collision, random, warnings, time, os
try: 
    import cPickle as pickle
except:
    import pickle as pickle

def make_cg_dna(dna_segnames, dna_resids, bp_per_bead, aa_all, frame=0):
    dna1 = dna_segnames[0]
    dna2 = dna_segnames[1]
    resid1 = dna_resids[0]
    resid2 = dna_resids[1]
    
    # check the input
    assert np.abs(resid1[1]-resid1[0]) == np.abs(resid2[1]-resid2[0]), (
        "number of paired bases in DNA strands are not the same")

    bps = np.abs(resid1[1]-resid1[0])+1
    nbeads = int(np.round(bps/float(bp_per_bead)))
    #print 'nbeads = ', nbeads

    # load the dna into its own sasmol object
    aa_dna = sasmol.SasMol(0)
    dna_filter = '((moltype[i] == "dna" or moltype[i] == "rna"))' # RNA problem
    error, aa_dna_mask = aa_all.get_subset_mask(dna_filter)
    error = aa_all.copy_molecule_using_mask(aa_dna, aa_dna_mask, frame) # select dna
    natomsD = aa_dna.natoms()
    print 'dna atoms = ', natomsD

    # populate the cg_dna sasmol object with semi-arbitrary atom properties
    cg_dna = sasmol.SasMol(0)
    basis_filter = ("((name[i] == 'N1') and (resid[i] <= "+str(nbeads)+") and "
                    "(segname[i] == '"+dna1+"'))")
    error, mask = aa_dna.get_subset_mask(basis_filter)   #; print mask
    s = np.sum(mask)
    assert s == nbeads, ("\n>>> ERROR!!! fail to create correct number of dummy"
                         " beads: expected %d but made %d\n") %(nbeads, s)
    error = aa_dna.copy_molecule_using_mask(cg_dna, mask, frame)

    # initialize bead coordinates and other variables
    cg_coor = np.zeros((1, nbeads, 3))
    r1a = resid1[0]
    r2a = resid2[0]
    all_beads = []       # list of the all atom sasmol object for each bead
    dna_bead_masks = []  # list of the all atom masks for each bead

    print "coarse-graining the DNA, this may take awhile..."
    tic = time.time()
    for j in xrange(nbeads):

        bead = sasmol.SasMol(0)
        if bp_per_bead > 1:
            if j+1 == nbeads:
                bp_per_bead = bps - bp_per_bead*j  # account for remainder

            # Get the atoms from DNA strand 1
            if resid1[0] < resid1[1]:
                r1b = r1a + bp_per_bead  # r1b > r1a
                basis_filter1 = ("((resid[i] > "+str(r1a-1)+" and "
                                 "resid[i] < "+str(r1b)+") and "
                                 "(segname[i]=='"+dna1+"')) or ")
            else:
                r1b = r1a - bp_per_bead  # r1b < r1a
                basis_filter1 = ("((resid[i] > "+str(r1b)+" and "
                                 "resid[i] < "+str(r1a+1)+") and "
                                 "(segname[i]=='"+dna1+"')) or ")
    
            # Get the atoms from DNA strand 2
            if resid2[0] < resid2[1]:
                r2b = r2a + bp_per_bead  # r2b > r2a
                basis_filter2 = ("((resid[i] > "+str(r2a-1)+" and "
                                 "resid[i] < "+str(r2b)+") and "
                                 "(segname[i]=='"+dna2+"'))")
            else:
                r2b = r2a - bp_per_bead   # r2b < r2a
                basis_filter2 = ("((resid[i] > "+str(r2b)+" and "
                                 "resid[i] < "+str(r2a+1)+") and "
                                 "(segname[i]=='"+dna2+"'))")
        else:
            # Get the atoms from DNA strand 1
            if resid1[0] < resid1[1]:
                r1b = r1a + bp_per_bead  # r1b > r1a
            else:
                r1b = r1a - bp_per_bead  # r1b < r1a
            
            basis_filter1 = ("((resid[i] == "+str(r1a)+") and "
                             "(segname[i]=='"+dna1+"')) or ")
    
            # Get the atoms from DNA strand 2
            if resid2[0] < resid2[1]:
                r2b = r2a + bp_per_bead  # r2b > r2a
            else:
                r2b = r2a - bp_per_bead   # r2b < r2a

            basis_filter2 = ("((resid[i]== "+str(r2a)+") and "
                             "(segname[i]=='"+dna2+"'))")
            
        basis_filter = basis_filter1+basis_filter2

        # create a mask to select the atoms for the bead
        error, mask = aa_dna.get_subset_mask(basis_filter)

        # store the mask for the reverse coarse-graining
        dna_bead_masks.append(mask)
        error = aa_dna.copy_molecule_using_mask(bead, mask, frame)

        com = bead.calccom(frame)
        cg_coor[0, j, :] = com  # a list of the com coordinate for each bead

        # calculate atomic coodinates using the bead com as the origin
        bead.center(0)
        all_beads.append(bead)
        
        # setup for next iteration
        r1a = r1b
        r2a = r2b
    
    cg_dna.setCoor(cg_coor)  #; print cg_dna._coor[frame, :, 0]

    vecXYZ = np.zeros((3, nbeads, 3))
    vecXYZ[0] = [1, 0, 0]
    vecXYZ[1] = [0, 1, 0]
    vecXYZ[2] = [0, 0, 1]
    toc = time.time() - tic
    print 'DNA coarse-graining took %f seconds' % toc
    return (aa_dna, aa_dna_mask, cg_dna, all_beads, dna_bead_masks, vecXYZ)

def make_cg_pro(aa_all, pro_groups, frame=0):
    # load the protein into its own sasmol object
    aa_pro = sasmol.SasMol(0)
    error, aa_pro_mask = aa_all.get_subset_mask('((moltype[i] == "protein"))')
    error = aa_all.copy_molecule_using_mask(aa_pro, aa_pro_mask, frame)
    print 'protein atoms = ', aa_pro.natoms()

    # coarse-grain the proteins
    print "coarse-graining the protein..."
    tic = time.time()
    cg_pro = sasmol.SasMol(0)
    error, aa2cg_pro_mask = aa_pro.get_subset_mask("(name[i] == 'CA')")
    error = aa_pro.copy_molecule_using_mask(cg_pro, aa2cg_pro_mask, frame)

    # create a sasmol object for each protein move group (pro_groups) and 
    # each group of all atom proteins
    cg_pgroup_masks = []  # masks to separate the cg_protein into NCP groups
    aa_pgroup_masks = []  # masks to separate the aa_protein into NCP groups
    all_proteins = []
    (pre, post, btwn) = ("((segname[i]=='", "'))", " or ")
    for group in pro_groups:
        aa_pro_group = sasmol.SasMol(0)
        if len(group) > 0:
            basis_filter = ""
            for seg in group:
                basis_filter += pre + seg + post + btwn
            basis_filter = basis_filter[0:-4]  # remove the last btwn
            error, cg_pgroup_mask = cg_pro.get_subset_mask(basis_filter)
            error, aa_pgroup_mask = aa_pro.get_subset_mask(basis_filter)
            error = aa_pro.copy_molecule_using_mask(aa_pro_group,
                                                    aa_pgroup_mask, frame)
        else:
            cg_pgroup_mask = np.zeros(sum(aa2cg_pro_mask), dtype='int32')
            aa_pgroup_mask = np.zeros(len(aa2cg_pro_mask), dtype='int32')

        cg_pgroup_masks.append(cg_pgroup_mask)
        aa_pgroup_masks.append(aa_pgroup_mask)
        all_proteins.append(aa_pro_group)

    # combined the protein groups into move groups
    move_masks = []
    for i in xrange(len(cg_pgroup_masks)):
        move_masks.append(np.copy(cg_pgroup_masks[i]))
    for i in xrange(len(move_masks)):
        for j in xrange(i+1, len(move_masks)):
            move_masks[i] += move_masks[j]
            
    toc = time.time() - tic
    print 'Protein coarse-graining took %f seconds' % toc
    
    return (aa_pro, aa_pro_mask, cg_pro, cg_pgroup_masks, aa_pgroup_masks, 
            all_proteins, move_masks)

def is_bead_flexible(flex_resids, nbeads, resid1, bp_per_bead):
    ## Flexible part
    # create an array (groupid) labeling the group for each flexible cg-DNA bead
    flexResids = [item for sublist in flex_resids for item in sublist]
    groupid = np.zeros(len(flexResids))
    link = np.zeros((nbeads-1, 2), dtype=int)  # links between DNA beads
    # if a bead is designated flexible both it's links are flexible

    n_start = 0
    for i in xrange(len(flex_resids)):
        n_end = len(flex_resids[i]) + n_start
        groupid[n_start:n_end] = np.ones(len(flex_resids[i])) * i
        n_start = n_end
    
    flexResids = np.concatenate(([flexResids], [groupid]),0).T    

    r1a = resid1[0]
    for j in xrange(nbeads):
        if resid1[0] < resid1[1]:
            r1b = r1a + bp_per_bead  # r1b > r1a
        else:
            r1b = r1a - bp_per_bead  # r1b < r1a            
        ## check residues from DNA strand 1 to see if this is a flexible bead
        for i in xrange(r1a, r1b):
            if i in flexResids[:, 0]:
                (r, c) = np.where(flexResids==i)
                # this stores the assigned group in the temp var: beadGroup
                beadGroup = flexResids[r[0], 1]
                if j==0:                         # first bead
                    link[j] = [1, beadGroup]
                elif j==nbeads-1:                # last bead
                    link[j-2] = [1, beadGroup]   # should this be j-1?  8/22/14
                else:  #all other beads, this is only 2 links in python
                    link[j-1:j+1] = [1, beadGroup]
                break  # bead defined as flexible, stop looping over residues        

        # setup for next iteration
        r1a = r1b

    # position the bead at its calculated com
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
    return beadgroups, trialbeads


def recover_aaPro_model(aa_pgroup_masks, cg_pgroup_masks, cg_pro, all_proteins,
                        aa_pro):
    """
    Recover the all-atom representation of coarse-grained proteins
    """
    for (i, mask) in enumerate(aa_pgroup_masks):
        if sum(mask) > 0:
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
            error, aa2cg_mask = all_proteins[i].get_subset_mask(
                "(name[i] == 'CA')")
            error = all_proteins[i].copy_molecule_using_mask(cg_pro_group_orig,
                                                             aa2cg_mask, 0)
            # get the sub_2 inputs for the align function
            com_sub_2 = cg_pro_group_orig.calccom(0)
            cg_pro_group_orig.center(0)
            coor_sub_2 = cg_pro_group_orig.coor()[0]

            all_proteins[i].align(0, coor_sub_2, com_sub_2,
                                  coor_sub_1, com_sub_1)
            error = aa_pro.set_coor_using_mask(all_proteins[i], 0, mask)

def recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, allBeads, masks):
    """
    Recover the all-atom representation of coarse-grained DNA
    """
    error =[]
    cg_natoms = np.copy(cg_dna.natoms())
    coor = np.copy(cg_dna.coor()[0])

    # split vecXYZ into three matrices
    vecX = np.copy(vecXYZ[0])
    vecY = np.copy(vecXYZ[1])
    vecZ = np.copy(vecXYZ[2])

    for i in xrange(cg_natoms):
        # R would align the rotated coordinates back to where they originated
        R = align2xyz(vecX[i, :], vecY[i, :], vecZ[i, :])

        # M will rotate cooridates from their original reference frame to the
        # rotated refrenece frame of the bead
        # then translate the com to the bead's com
        M = R.transpose()
        M[3, :3] = coor[i, :] # insert the com translation into the M matrix

        r, c = allBeads[i].coor()[0].shape
        beadCoor = np.ones((r, 4))           # initialize the beadCoor matrix
        origCoor = np.ones((1, r, 3))
        origCoor = allBeads[i].coor()
        beadCoor[:, :3] = np.copy(origCoor[0])
        beadCoor = np.dot(beadCoor, M)

        tmpCoor = np.zeros((1, r, 3))        # initialize the tmpCoor
        tmpCoor[0] = beadCoor[:, :3]
        allBeads[i].setCoor(tmpCoor)

        #recombine all the beads back into one pdb
        e = aa_dna.set_coor_using_mask(allBeads[i], 0, masks[i])
        error.append(e)

        #reset the bead coordinates for the next iteration
        allBeads[i].setCoor(origCoor)

    return error

def align2xyz(vecX, vecY, vecZ):

    tmp_coor = np.zeros((2, 4))
    tmp_coor[:, 3] = 1
    tmp_coor[1, 0:3] = vecZ
    A1 = align2z(tmp_coor)

    newX = np.dot(vecX, A1[0:3, 0:3]) #; print 'newX = ', newX
    #newY = np.dot(vecY, A1[0:3, 0:3]) ;# print 'newY = ', newY
    assert newX[2] < 10e-5, ("ERROR!!! z-component of newX is not zero and it "
                             "should be")

    thetaZ_x = -np.arctan2(newX[1], newX[0])
    #thetaZ_y = -np.arctan2(newY[1], -newY[0])

    A2 = rotate4x4('z', thetaZ_x)  #; print 'A2 = ', A2

    A = np.dot(A1, A2)

    newY = np.dot(vecY, A[0:3, 0:3]) #; print 'finalY = ', newY
    assert newY[0]+newY[2] < 1+10e-5, ("ERROR!!! newY is not aligned to the "
                                       "y-axis and it should be")

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
    align the axis connecting the first 2 coordinates of a 1x4 array of
    coodinate vectors to the z-axis
    '''
    A = np.eye(4, dtype=np.float)
    #s Axz  = np.eye(4, dtype=np.float)
    #s Az   = np.eye(4, dtype=np.float)
    assert all(coor4[0] == [0., 0., 0., 1., ]), ("coordinates passed to align2z"
                    "were not translated to the origin")

    if coor4.shape > (1, 4):
        (u, v, w) = coor4[1, 0:3]
        small = 1E-14 # small limit to make sure not dividing by zero

        d1 = np.sqrt(u**2+v**2)
        d2 = np.sqrt(u**2+v**2+w**2)
        if d1 > small:
            (A[0][0], A[0][1], A[0][2]) = (u/d1*w/d2, -v/d1, u/d2)
            (A[1][0], A[1][1], A[1][2]) = (v/d1*w/d2, u/d1, v/d2)
            (A[2][0], A[2][1], A[2][2]) = (-d1/d2, 0, w/d2)

        #s # align to the x-z plane
        #s if v > small:           # check if already aligned to xz-plane
        #s         (Axz[0][0], Axz[0][1]) = (u/d1, -v/d1)
        #s         (Axz[1][0], Axz[1][1]) = (v/d1, u/d1)
        #s
        #s # align to the z-axis
        #s if d1 > small:
        #s         (Az[0][0], Az[0][2]) = (w/d2, d1/d2)
        #s         (Az[2][0], Az[2][2]) = (-d1/d2, w/d2)
    else:
        print 'no point to align'
        #s A = np.dot(Axz, Az) #;   print 'A3=\n', A
    return A

def beadRotate(coor3, vecXYZ, thetas, n_soft, p_coor3):
    '''
    this function is designed to generate a modified version of the input
    coordinates (coor3)
    1) translate all coordinates so the first one is at the origin
    2) orient the second coordinate along the z-axis
    3) performs rotations using the input thetas averaged over n_soft beads
    4) reverse the orientation to the z-axis
    5) reverse transaltion of all coordinates so the first is where it started
    '''
    (natoms, col) = coor3.shape

    # initialize vector arrays for coordinates and orientation vectors
    # changing them from 3 component vectors into 4 component vectors to
    # incorporate transaltions into the matrix math
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
    # This is intended to be multiplied from the right (unlike standard matrix
    # multiplication) so as not to require transposing the coordinate vectors.
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
    (R[2][0], R[2][1], R[2][2]) = ( sx*sz+cx*sy*cz, cx*sy*sz-sx*cz, cx*cy )

    # verbose method of defining the rotation matrices
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

    (T0, Ti0) = move2origin(coor4)  # Ti0 != T0, negative off diag elements
    coor4 = np.dot(coor4, T0)  #; print 'moved to origin coor:\n', coor4
    A = align2z(coor4)
    AR = np.dot(A, R)

    # first of n_soft rotations
    coor4 = np.dot(coor4, AR)

    # the coarse grained beads local coordinates should not be translated,
    # only rotated
    X = np.dot(X, AR)
    Y = np.dot(Y, AR)
    Z = np.dot(Z, AR)

    p_coor4 = np.dot(np.dot(p_coor4, T0), AR)

    # remaining of n_soft rotations
    for i in xrange(1, n_soft):
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

        # coarse grained beads local coordinates should not be translated,
        # only rotated
        X[i:] = np.dot(X[i:], R)
        Y[i:] = np.dot(Y[i:], R)
        Z[i:] = np.dot(Z[i:], R)


    coor4 = np.dot(coor4, A.transpose())
    coor4 = np.dot(coor4, Ti0)
    X = np.dot(X, A.transpose())
    Y = np.dot(Y, A.transpose())
    Z = np.dot(Z, A.transpose())

    p_coor4 = np.dot(p_coor4, A.transpose())
    p_coor4 = np.dot(p_coor4, Ti0)

    vecXYZ[0, 1:] = X[1:, 0:3]
    vecXYZ[1, 1:] = Y[1:, 0:3]
    vecXYZ[2, 1:] = Z[1:, 0:3]


    # this returns the modified positions and orientations for all but the
    # first bead or reference bead
    return (coor4[1:, 0:3], vecXYZ[:, 1:], p_coor4[:, 0:3])

def checkU(coor):
    '''
    module to generate and check the u-vectors pointing from bead to bead
    (used to make sure they were are all the same length but this did not work
    with real coarse-grained DNA)
    returns the u-vectors and their average magnitude (l)
    '''
    erLim = 1e-3
    (r, c) = coor.shape
    u = coor[1:, :]-coor[:r-1, :]                  # u-vectors between beads
    lu = np.sqrt(u[:, 0]**2+u[:, 1]**2+u[:, 2]**2) # magnitues of u-vectors
    l = np.mean(lu)                                # average distance btwn bead

    #print '\n\ncoor:\n', coor
    #print 'u-vectors:\n', u
    #print 'magnitude of u-vectors:\n', lu
    #print 'average magnitude of u =', l
    #print '\n\n'

    # Commented this out because it did not work for real DNA because base pairs
    # did not always start evenly spaced apart
    #s test = np.abs(lu-l) > erLim  # check if any of the l-lengths are different
    #s if test.any():
    #s         print 'ERROR: the beads are not uniformly spaced'
    #s         print 'u = \n', u
    #s         print test
    #s         # print 'difference = ', u[:, 2] - l
    #s         print 'distance = ', lu

    return (u, l)

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

def p_energy_wca(w, coor, wca0, trial_bead):
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
            rij = np.sqrt((ri[0]-rj[0])**2. +
                          (ri[1]-rj[1])**2. +
                          (ri[2]-rj[2])**2.)
            #s print '(2^(1/6)*w, rij) = ', (test, rij)
            if rij < test:
                wca1[i, j] = (w/rij)**12.-(w/rij)**6.+0.25
            else:
                wca1[i, j] = 0
    res = 4*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def f_energy_wca(w, coor, wca0, trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN
    '''

    wca1 = np.copy(wca0)
    # calculate the distance between moved DNA beads and all other DNA beads
    # increment trial bead by one because it will index from 1 instead of 0
    wca1 = collision.wca_d(coor, trial_bead+1, w, wca1)

    res = 4.*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def f_overlap1(coor, cutoff):
    '''
    this function checks for overlap using FORTRAN
    '''
    # calculate the distance between moved DNA beads and all other DNA beads
    check = collision.overlap1(coor, cutoff)

    return check

def f_overlap2(coor_a, coor_b, cutoff):
    '''
    this function checks for overlap using FORTRAN
    '''
    if len(coor_a) > 0 and  len(coor_b) > 0:
        # calculate the distance between moved DNA beads and all other DNA beads
        check = collision.overlap2(coor_a, coor_b, cutoff)
    else:
        check = 0

    return check

def mask2ind(mask):
    '''
    convert a mask to an array of indicies
    '''

    mask_indices = np.nonzero(mask*np.arange(1, len(mask)+1))[0]

    return mask_indices

def dna_mc(ARGS, cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, lp, trialbeads,
           beadgroups, group_masks, all_beads, dna_bead_masks, aa_pgroup_masks,
           cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask,
           dna_type='b'):
    '''
    this function perform nsteps Monte-Carlo moves on the cg_dna
    '''

    timestr = time.strftime("%y%m%d_%H%M%S_") # prefix for output files
    all_dcd_name = timestr + ARGS.pdb[:-4] + '.dcd'
    aa_all_dcd_out = aa_all.open_dcd_write(all_dcd_name)
    
    if False:
        aa_all.send_coordinates_to_vmd(2222,0)
    
    # create the coarse-grained DNA and protein dcd and pdb files
    cg_dna_dcd_name = timestr + 'cg_dna.dcd'
    cg_pro_dcd_name = timestr + 'cg_pro.dcd'    
    cg_dna_dcd_out = cg_dna.open_dcd_write(cg_dna_dcd_name)
    cg_pro_dcd_out = cg_pro.open_dcd_write(cg_pro_dcd_name)    
    cg_dna.write_dcd_step(cg_dna_dcd_out, 0, 1)
    cg_pro.write_dcd_step(cg_pro_dcd_out, 0, 1)    
    cg_dna.write_pdb(timestr + 'cg_dna.pdb', 0, 'w')
    cg_pro.write_pdb(timestr + 'cg_pro.pdb', 0, 'w')    
    
    # create a dummy sasmol object for the 3 orientation vectors for each bead
    # will write these out to dcd files to store the coordinates along the way
    vecX_mol = sasmol.SasMol(0)
    vecY_mol = sasmol.SasMol(0)
    vecZ_mol = sasmol.SasMol(0)    
    error, mask = cg_dna.get_subset_mask('(all)')
    error = cg_dna.copy_molecule_using_mask(vecX_mol,mask,0)
    error = cg_dna.copy_molecule_using_mask(vecY_mol,mask,0)
    error = cg_dna.copy_molecule_using_mask(vecZ_mol,mask,0)
    vecX_mol.setCoor(np.array([vecXYZ[0]])) # the np.array recast these so they 
    vecY_mol.setCoor(np.array([vecXYZ[1]])) # do not update with vecXYZ 
    vecZ_mol.setCoor(np.array([vecXYZ[2]]))
    vecX_dcd_name = timestr + 'vecX.dcd'
    vecY_dcd_name = timestr + 'vecY.dcd'
    vecZ_dcd_name = timestr + 'vecZ.dcd'
    vecX_dcd_out = vecX_mol.open_dcd_write(vecX_dcd_name)
    vecY_dcd_out = vecY_mol.open_dcd_write(vecY_dcd_name)    
    vecZ_dcd_out = vecZ_mol.open_dcd_write(vecZ_dcd_name)        
    vecX_mol.write_dcd_step(vecX_dcd_out, 0, 1)
    vecY_mol.write_dcd_step(vecY_dcd_out, 0, 1)   
    vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, 1)   
    
    # initialize variables for each run
    steps_from_0 = np.zeros(ARGS.nsteps, dtype='int64')
    xyz = np.copy(vecXYZ)
    d_coor = np.copy(cg_dna.coor()[0]) # unique memory for each
    p_coor = np.copy(cg_pro.coor()[0]) # unique memory for each

    (u, l) = checkU(d_coor) # vectors between beads u, and average distance l
    #s print "(u, l) =", (u, l) # debug info
    lpl = lp/l  # setup the presistence length paramater

    dna_energy_width = {'a': 0, 'b': 46., 'z': 0} # yet to use a, and z type dna
    w = dna_energy_width[dna_type.lower()]
    if w > l:
        w = np.floor(l)
        # print '~~~ %.2f > %.2f ~~~~~~~~~~~~~~~~~~~~~~~~' % (w, l)
        print ('>>> setting chain width (w) to %d (chain width < distance' %w,
               ' btwn beads)')

    dna_diam = {'a': 25.5, 'b': 23.7, 'z': 18.4}
    dna_bead_radius = 4.5

    pro_bead_radius = 1.0 # 2A min seperation of CA atoms in database

    pro_pro_test = pro_bead_radius + pro_bead_radius
    dna_pro_test = dna_bead_radius + pro_bead_radius

    # calculate the energy of the starting positions
    wca0 = np.zeros((cg_dna.natoms(),cg_dna.natoms()))
    Ub0 = energyBend(lpl, u, l)

    if ARGS.f_collide:
        (Uwca0, wca0) = f_energy_wca(w, d_coor, wca0, 0)
    else:
        (Uwca0, wca0) = p_energy_wca(w, d_coor, wca0, 0)

    U_T0 = Ub0 + Uwca0
    # print '(Ub0, Uwca0, Ub0/U_T0, Uwca0/U_T0) = ', (Ub0, Uwca0, Ub0/U_T0, 
    #                                                 Uwca0/U_T0)
    n_accept   = 0 # total times configuration was accepted
    n_reject   = 0 # total times configuration was rejected
    n_written  = 0 # total times dcd write has been called
    fail_tally = 0 # number of times failed for particular iteration
    n_from_reload = 0 # number of stps since last reload
    n_reload = [0]  # listt containing the i_goback values
    
    # this should not actually be >=, come back to this
    assert np.size(ARGS.theta_max) - 1 >= np.max(beadgroups), (
        'each group needs its own theta_max: %d < %d'
        % (np.size(ARGS.theta_max) - 1, np.max(beadgroups) ))

    # Main MC loop #
    while n_accept < ARGS.nsteps:

        # Choose a bead to rotate
        trial_bead = trialbeads[int((trialbeads.size)*np.random.random())]        

        # Determine rotation to perform
        theta_max = ARGS.theta_max[beadgroups[trial_bead]]
        thetaZ_max = 0 * np.float(theta_max) # option to scale thetaZ separatly
        #thetaZ_max = np.float(theta_max) # option to scale thetaZ separatly
        thetaZ = 2 * thetaZ_max * np.random.random() - thetaZ_max
        thetaX = 2 * theta_max  * np.random.random() - theta_max
        thetaY = 2 * theta_max  * np.random.random() - theta_max
        thetaXYZ = [thetaX/ARGS.n_soft, thetaY/ARGS.n_soft, thetaZ/ARGS.n_soft]
        # print theta_max, thetaXYZ

        if  len(group_masks) == 0 or beadgroups[trial_bead] == len(group_masks):
            # Only DNA will be moving, create place-holder dummy coordinates
            p_coor_rot = np.zeros((0, 3))
        else:
            p_mask = group_masks[beadgroups[trial_bead]]
            p_ind_rot = mask2ind(p_mask)
            p_ind_fix = mask2ind(-(p_mask-1))
            p_coor_rot = p_coor[p_ind_rot]
            p_coor_fix = p_coor[p_ind_fix]

        (d_coor[trial_bead:], xyz[:, trial_bead:], p_coor_rot) = beadRotate(
            d_coor[trial_bead-1:], xyz[:, trial_bead-1:], thetaXYZ, ARGS.n_soft,
            p_coor_rot) # generate a newly rotated model

        # store the rotated protein coordinates
        if beadgroups[trial_bead] < len(group_masks):
            p_coor[p_ind_rot] = p_coor_rot

        # calculate the change in energy (dU) and the boltzman factor (p)
        (u, l) = checkU(d_coor)
        Ub1 = energyBend(lpl, u, l)

        # ~~~~ DNA interaction energy  ~~~~~~#
        if ARGS.f_collide:
            (Uwca1, wca1) = f_energy_wca(w, d_coor, wca0, trial_bead)
        else:
            (Uwca1, wca1) = p_energy_wca(w, d_coor, wca0, trial_bead)
            print "python wca calculator deprecated"

        U_T1 =  Ub1 + Uwca1
        dU = U_T1 - U_T0

        with warnings.catch_warnings():
            warnings.filterwarnings('error') # need this for np warnings
            try:
                p = np.exp(-dU)
                # print '\n(Ub1, Uwca1) = ', (Ub1, Uwca1) 
                # print '(Ub1/U_T1, Uwca1/U_T1) = ', (Ub1/U_T1, Uwca1/U_T1)
                # print '(p, dU) = ', (p, dU)
            except Warning:
                if dU > 99:
                    p =  0
                    #s print 'energy was large, setting probability to 0'
                elif dU < 0:
                    p =  1
                    #s print 'energy was negative, setting probability to 1'
                else:
                    print 'Warning: ~~> unclear OverflowError <~~ dU = ', dU
                    print 'not sure where the error originated from'

        test = np.random.random()
        collision = 0

        if test >= p:
            dna_pass = False
            # print 'step failed because of DNA energy'
        else:
            dna_pass = True

            # now check for collisions
            if len(p_coor_rot) > 0:   # only if proteins were rotated
                # ~~~~ Check for overlap, DNA-protein or protein-protein ~~~~~~#
                d_coor_fix = d_coor[trial_bead:]
                d_coor_rot = d_coor[:trial_bead]

                # check for protein-protein overlap
                if 1 == f_overlap2(p_coor_rot, p_coor_fix, pro_pro_test):
                    print 'Protein-Protein'
                    #print 'collision, set p=0'
                    collision = 1
                
                print 'currently ignoring DNA-protein overlap'
                # # check for DNA-protein overlap
                # elif 1 == f_overlap2(p_coor_rot, d_coor_fix, dna_pro_test):
                    # print 'Potein-DNA (rot-fix)'
                    # #print 'collision, set p=0'
                    # collision = 1
                    # print 'ignoring this for now'

                # elif 1 == f_overlap2(p_coor_fix, d_coor_rot, dna_pro_test):
                    # print 'Potein-DNA (fix-rot)'
                    # #print 'collision, set p=0'
                    # collision = 1

                if collision == 1:
                    print 'failed because of collision'

        if dna_pass and collision == 0:
            n_from_reload += 1
            steps_from_0[n_accept] = n_from_reload + n_reload[-1]
            n_accept += 1                      # increment accept counter
            # cg_dna.setCoor(d_coor) # <-- DO NOT use setCoor, want uniuqe mem
            # cg_pro.setCoor(p_coor) # <-- DO NOT use setCoor, want uniuqe mem            
            cg_dna.setCoor(np.array([(d_coor)])) # update dna coordinates
            cg_pro.setCoor(np.array([(p_coor)])) # update protein coordinates
            vecXYZ = np.copy(xyz)              # update dna orientations
            vecX_mol.setCoor(np.array([vecXYZ[0]])) # independent of vecXYZ[0]
            vecY_mol.setCoor(np.array([vecXYZ[1]])) # independent of vecXYZ[1]
            vecZ_mol.setCoor(np.array([vecXYZ[2]])) # independent of vecXYZ[2]
            
            wca0 = np.copy(wca1)               # update DNA WCA energy        
            U_T0 = U_T1                        # update total energy

            #print output regarding trial
            print "trial_bead(%3d) = %2d\t failed attempts = %2d" % (n_accept, 
                                                        trial_bead, fail_tally)
            fail_tally = 0                     # reset fail_tally

            # write out the accepted configuration for go-back use
            if ARGS.goback > 0:
                # these are incremented by one because the original coordinates 
                # are saved (that is not the case for aa_all)
                cg_dna.write_dcd_step(cg_dna_dcd_out, 0, n_written+1)
                cg_pro.write_dcd_step(cg_pro_dcd_out, 0, n_written+1)                
                vecX_mol.write_dcd_step(vecX_dcd_out, 0, n_written+1)
                vecY_mol.write_dcd_step(vecY_dcd_out, 0, n_written+1)   
                vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, n_written+1)   
            
            # recover an all atom representation and save coordinates to a dcd
            # this requires re-inserting the aa-coordinates which takes added 
            # time so only do when designated
            if 0 == n_accept % ARGS.n_dcd_write:
                # ~~recover aa-DNA~~
                error = recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, all_beads,
                                            dna_bead_masks)
                # ~~recover aa-Protein~~
                recover_aaPro_model(aa_pgroup_masks, cg_pgroup_masks, cg_pro,
                                    all_proteins, aa_pro)
                # ~~Combine aa Complete Structure~~
                aa_all.set_coor_using_mask(aa_pro, 0, aa_pro_mask)
                aa_all.set_coor_using_mask(aa_dna, 0, aa_dna_mask)
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written)
                
        else :
            # by default ARGS.goback is -1 so this returns FALSE without input
            if fail_tally == ARGS.goback:  
                i_goback = rewind(ARGS, n_accept, cg_dna_dcd_name,
                            cg_dna, cg_pro_dcd_name, cg_pro, vecX_dcd_name, 
                            vecX_mol, vecY_mol, vecY_dcd_name, vecZ_mol, 
                            vecZ_dcd_name, vecXYZ)

                d_coor = np.copy(cg_dna.coor()[0]) # reset the dna coordinates

                # reset the reference energy
                (u, l) = checkU(d_coor) 
                (Uwca0, wca0) = f_energy_wca(w, d_coor, wca0, 0)
                Ub0 = energyBend(lpl, u, l)
                U_T0 =  Ub0 + Uwca0
                
                n_from_reload = 0
                n_reload.append(steps_from_0[i_goback-1])
                fail_tally = 0 # reset the fail counter
            else:
                fail_tally += 1                 # increment bead reject counter 
                n_reject += 1                   # increment total reject counter
                d_coor = np.copy(cg_dna.coor()[0]) # reset the dna coordinates
                
            p_coor = np.copy(cg_pro.coor()[0]) # reset the protein coordinates
            xyz = np.copy(vecXYZ)              # reset the dna orientations

            # save previous coordinates again
            if not ARGS.keep_unique:
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written )
                cg_dna.write_dcd_step(cg_dna_all_dcd_out, 0, n_written + 1)
                cg_pro.write_dcd_step(cg_pro_all_dcd_out, 0, n_written + 1)                

    aa_all.close_dcd_write(aa_all_dcd_out)
    
    cg_dna.close_dcd_write(cg_dna_dcd_out) #uncomment if wanting to keep
    cg_pro.close_dcd_write(cg_pro_dcd_out) #uncomment if wanting to keep
    # os.remove(timestr + 'cg_dna.pdb')
    # os.remove(timestr + 'cg_pro.pdb')    
    # os.remove(cg_dna_dcd_name) # remove/comment to keep the cg dna coor
    # os.remove(cg_pro_dcd_name) # remove/comment to keep the cg pro coor   
    os.remove(vecX_dcd_name) # remove/comment to keep the cg pro coor   
    os.remove(vecY_dcd_name) # remove/comment to keep the cg pro coor   
    os.remove(vecZ_dcd_name) # remove/comment to keep the cg pro coor       

    if ARGS.goback > 0:
        np.savetxt(timestr+'n_from_0.txt', steps_from_0, fmt='%d')
        
    print "accepted %d moves" % n_accept
    print "rejected %d moves" % n_reject

    # print n_reload
    # print steps_from_0

def rewind(ARGS, n_accept, cg_dna_dcd_name, cg_dna, cg_pro_dcd_name, cg_pro, 
           vecX_dcd_name, vecX_mol, vecY_mol, vecY_dcd_name, vecZ_mol,
           vecZ_dcd_name, vecXYZ):
    
    if ARGS.seed > 0:
        np.random.seed(ARGS.seed)
        i_goback = int(n_accept*np.random.random())
    else:
        i_goback = int(n_accept*np.random.random())
        
    # i_goback = 0 represents the original structure this is the 1st frame of 
    # the cg_dna_dcd and cg_pro_dcd so increment by 1
        
    cg_dna.read_single_dcd_step(cg_dna_dcd_name, i_goback + 1)
    cg_pro.read_single_dcd_step(cg_pro_dcd_name, i_goback + 1) 
    vecX_mol.read_single_dcd_step(vecX_dcd_name, i_goback + 1)
    vecY_mol.read_single_dcd_step(vecY_dcd_name, i_goback + 1)   
    vecZ_mol.read_single_dcd_step(vecZ_dcd_name, i_goback + 1)                                     
    vecXYZ[0] = vecX_mol.coor()[0]
    vecXYZ[1] = vecY_mol.coor()[0]
    vecXYZ[2] = vecZ_mol.coor()[0]                

    if i_goback == 0:
        print '\n~~~ reloaded original coordinates ~~~\n'
    else:	
        print '\n~~~ reloaded accepted coordinates #%d ~~~\n' % i_goback
        
    return i_goback

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

    longCoor = np.zeros((1, natoms, 3), np.float) # initialize long DNA position
    longCoor[0][:, 2] = range(natoms) # set the z-values to the array index
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

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'Generate modified DNA or DNA-protein structures'
        #epilog = 'no epilog found'
    )
    parser.add_argument("--seed", default=0, type=int,
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
    parser.add_argument("-tm", "--theta_max", nargs='+', type=np.float,
        help = ("Max rotation angle for the x and y rotation; "
                "max z scales with this paramater (curently 1x but code is "
                "designed to allow max z to scale differently from max x and "
                "y); default = 5"), default=[25])
    parser.add_argument("-ns", "--n_soft", default=1, type=int,
        help = ("Number of bead over which to apply any rotation (i.e. soften "
                "the bending); default = 1") )
    parser.add_argument("-bp", "--bp_per_bead", default=1, type=int,
        help = "Number of bp for each coarse-grained bead; default = 1")
    parser.add_argument("-fl", "--flex", default=False, action="store_true",
        help = ("Generate output file with the flexible resids formatted into "
                "2 columns.  The first row contains the segnames for the 2 "
                "DNA strands; default = False"))

    parser.add_argument("--size", default='med', type=str,
                        help="selection for new_c11h5.pdb flexible residues")
    parser.add_argument("--regions", default='1-5', type=str,
                        help="selection for new_c11h5.pdb flexible regions")

    parser.add_argument("-rp", "--rm_pkl", default=False, action="store_true",
                        help="remove the pkl file containing cg parameters")

    parser.add_argument("-f", dest="f_collide", action='store_true',
        help="use fortran module for calculating collisions; the default")
    parser.add_argument("-py", dest="f_collide", action='store_false',
        help="use python module for calculating collisions")
    parser.set_defaults(f_collide=True)

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

def write_flex_resids(all_beads, flex_beads, ARGS):
    '''
    Write the segnames names and resids of a list of sasmol objects to file
    Designed to be used for a list of sasmol objects representing coarse-grained
    DNA beads.
    These files can be read using the function 'read_flex_resids'
    '''
    dna1_bead_resids = []
    dna2_bead_resids = []
    dna_bead_resids = [dna1_bead_resids, dna2_bead_resids]
    dna1_bead_labels = []
    dna2_bead_labels = []
    dna_bead_labels = [dna1_bead_labels, dna2_bead_labels]

    bases_per_bead = np.zeros((len(flex_beads)), dtype='int')
    fewer_bp = []
    for (i, flex_bead) in enumerate(flex_beads):
        bead = all_beads[flex_bead]
        for (j, label) in enumerate(bead.segnames()):
            basis_filter = "((segname[i]=='"+label+"'))"
            error, dna_strand_mask = bead.get_subset_mask(basis_filter)
            dna = sasmol.SasMol(0)
            error = bead.copy_molecule_using_mask(dna, dna_strand_mask, 0)
            # dna1_bead_labels.append(all_beads[i].segnames()[0])
            dna_bead_labels[j].append(dna.segnames()[0])
            dna_bead_resids[j].append(dna.resids())
            bases_per_bead[i] += len(dna.resids())
        if bases_per_bead[i] < bases_per_bead[0]:
            fewer_bp.append(i)
            n_fewer = bases_per_bead[0]/2 - bases_per_bead[i]/2
            print '%d fewer bases in bead %d than the others' % (n_fewer, i + 1)
            print 'appending 0 to the flex file'
            for j in  xrange(n_fewer):
                dna1_bead_resids[i].append(0)
                dna2_bead_resids[i].append(0)

    dna1_label = list(set(dna1_bead_labels))
    dna2_label = list(set(dna2_bead_labels))
    assert 1 == len(dna1_label) & 1 == len(dna2_label), ("multiple segnames "
                                                         "within DNA strand")

    flex_id_out = open(ARGS.pdb[:-3]+'flex', 'w')
    flex_id_out.write('%s %s\n' % (dna1_label[0], dna2_label[0]))

    n_flex_bp = len(flex_beads)*bases_per_bead[0] / 2
    dna1_resids = np.reshape(np.array(dna1_bead_resids), n_flex_bp)
    dna2_resids_temp = np.reshape(np.array(dna2_bead_resids), n_flex_bp)
    dna2_resids_temp.sort()
    dna2_resids = dna2_resids_temp[::-1]
    for i in xrange(n_flex_bp):
        flex_id_out.write('%d %d\n' % (dna1_resids[i], dna2_resids[i]))
    flex_id_out.close()

def read_flex_resids(flex_file):
    '''
    Read flexible DNA resids from the file created using 'write_flex_resids'
    First line of file should have the DNA chain names
    Following lines should have the resids from each chain that are designated
    as flexible, for example:
    A B
    1 4
    2 3
    3 2
    4 1
    '''
    lines = [line.strip() for line in open(flex_file)]
    segnames = []
    segnames.append(lines[0].split()[0])
    segnames.append(lines[0].split()[1])
    flex_resids = np.genfromtxt(flex_file, dtype='int',delimiter=" ")[1:]
    return (segnames, flex_resids)

def main():

    # set the DNA properties parameters:
    lp = 530.     # persistence length  (lp = 530A)
    # defined in dna_mc
    #w = 46        # effective measure of dsDNA chain width in A (w = 46A)
    #DOI 10.1021/ma201277e used w=46, lp=530

    if ARGS.f_collide:
        print 'using fortran collision calculation'
    else:
        print 'using python collision calculation'

    if ARGS.keep_unique:
        print 'only keeping unique structures'
    else:
        print 'keeping all accepted structures'

    if ARGS.Llp:
        print 'simulating long DNA'
        Llp = ARGS.Llp    # Llp is the length in units of lp: Llp = L/lp
        L = Llp*lp  #
        (cg_dna, vecXYZ) = makeLongDNA(Llp) # use this to make long cgDNA
        all_atom_pdb = '%d_llp' % Llp

    elif ARGS.pdb:
        print 'loading pdb'
        all_atom_pdb = ARGS.pdb
        dna_resids =  []
        pro_groups =  []

        if  ARGS.pdb ==  'new_c11_tetramer.pdb':
            '''The C11 nucleosome tetramer with all the protein tails'''
            dna_segnames = ['DNA1', 'DNA2']
            dna_resids.append([1, 693])
            dna_resids.append([693, 1])
            # continuous flexible residues on the first DNA strand
            # recall that range(a, b) excludes upper lim: [a, b)
            flex_resids = [range(1, 31), range(167, 198), range(334, 365),
                 range(501, 532), range(667, 694)]
            pro_groups.append(['A0', 'B0', 'C0', 'D0',
                               'E0', 'F0', 'G0', 'H0'])
            pro_groups.append(['A1', 'B1', 'C1', 'D1',
                               'E1', 'F1', 'G1', 'H1'])
            pro_groups.append(['M1', 'N1', 'O1', 'P1',
                               'Q1', 'R1', 'S1', 'T1'])
            pro_groups.append(['M0', 'N0', 'O0', 'P0',
                               'Q0', 'R0', 'S0', 'T0'])
            ARGS.theta_max = [10, 10, 10, 10, 10]

        elif  ARGS.pdb ==  'new_c11h5.pdb':
            '''The C11 nucleosome tetramer with the gH5 linker'''
            dna_segnames = ['DNA1', 'DNA2']
            dna_resids.append([1, 693])
            dna_resids.append([693, 1])
            # continuous flexible residues on the first DNA strand
            #s flex_resids = [range(1, 31), range(167, 198), range(334, 365),
            #s     range(501, 532), range(667, 694)]
            small_range = [range(1, 10), range(177, 180), range(344, 352),
                 range(513, 517), range(680, 694)]
            inter_range = [range(1, 20), range(172, 189), range(339, 357),
                 range(507, 524), range(674, 694)]
            large_range = [range(1, 31), range(167, 198), range(334, 365),
                 range(501, 532), range(667, 694)]
            xlarge_range = [range(1, 36), range(162, 203), range(329, 370),
                 range(496, 537), range(662, 699)]

            if ARGS.size == 'large':
                flex_resids = large_range
            elif ARGS.size == 'med':
                flex_resids = inter_range
            elif ARGS.size == 'small':
                flex_resids = small_range
            elif ARGS.size == 'x-lrg':
                flex_resids = xlarge_range
            else:
                print 'unknown selection for size: ', size
                
            if ARGS.regions == '1-5':         # 1-5 option
                pass                          # flex_resids = flex_resids
            elif ARGS.regions == '234':       # 2-4 option
                flex_resids[0] = flex_resids[4] = []
            elif ARGS.regions == '2_4':       # 2, 4 option
                flex_resids[0] = flex_resids[4] = flex_resids[2] = []
            elif ARGS.regions == '23_':       # 2-3 option
                flex_resids[0] = flex_resids[4] = flex_resids[3] = []
            elif ARGS.regions == '2__':       # 2 option
                flex_resids[0] = flex_resids[4] = flex_resids[3] = flex_resids[2] = []
            elif ARGS.regions == '_3_':       # 3 option
                flex_resids[0] = flex_resids[4] = flex_resids[3] = flex_resids[1] = []
            else:
                print 'unknown selection for flex: ', ARGS.regions
                
            # manually selecting flex-resids
            #s flex_resids = [range(18,33), range(182,199), [], range(496,508), range(668,677)]
            flex_resids = [[], [], [], [], range(668,677)]
            

            print ('using regions ' + ARGS.regions + ' with a ' + ARGS.size + 
                   ' size -> flex_resids = \n', flex_resids)

            pro_groups.append(['A0', 'B0', 'C0', 'D0',
                               'E0', 'F0', 'G0', 'H0', 'H5N1'])
            pro_groups.append(['A1', 'B1', 'C1', 'D1',
                               'E1', 'F1', 'G1', 'H1', 'H5N2'])
            pro_groups.append(['M1', 'N1', 'O1', 'P1',
                               'Q1', 'R1', 'S1', 'T1', 'H5N3'])
            pro_groups.append(['M0', 'N0', 'O0', 'P0',
                               'Q0', 'R0', 'S0', 'T0', 'H5N4'])
            ARGS.theta_max = [2, 5, 2, 5, 2]

        elif ARGS.pdb == 'new_dsDNA60.pdb':
            # linker dna file
            dna_segnames = ['DNA1', 'DNA2']
            dna_resids.append([1, 60]) # DNA base pairing
            dna_resids.append([120, 61]) # DNA base pairing
            #flex_resids = [range(16, 45)]
            flex_resids = [range(1,15), range(46,60)]
            ARGS.theta_max = [5, 5, 5, 5, 5]
        elif ARGS.pdb == 'new_dsDNA.pdb':
            # linker dna file
            dna_segnames = ['DNA1', 'DNA2']
            dna_resids.append([1, 30]) # DNA base pairing
            dna_resids.append([30, 1]) # DNA base pairing
            flex_resids = [range(1, 31)]
        else:
            print "\n~~~ ERROR, unknow pdb file input ~~~\n"

        tic = time.time()

        (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, trialbeads, beadgroups, 
        move_masks, all_beads, dna_bead_masks, aa_pgroup_masks, cg_pgroup_masks, 
        all_proteins, aa_all, aa_pro_mask, aa_dna_mask) = get_cg_parameters(
            ARGS, flex_resids, pro_groups, dna_resids, dna_segnames)
        
        toc = time.time() - tic 
        print 'Total coarse-grain time = %f seconds' % toc        
        
    tic = time.time()     
    dna_mc(ARGS, cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, lp, trialbeads,
           beadgroups, move_masks, all_beads, dna_bead_masks, aa_pgroup_masks,
           cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask)
    toc = time.time() - tic; print 'run time =', toc, 'seconds'

    if None == ARGS.Llp and ARGS.flex:
        write_flex_resids(all_beads, trialbeads, ARGS)

def get_cg_parameters(ARGS, flex_resids, pro_groups, dna_resids, 
                      dna_segnames):
    '''
    This module is designed to generate the coarse-grained run parameters then
    save them to a pickle file for future use.  If the pickle file already 
    exists, it will load that file then compare the new input to see if all 
    parameters meet the current input.  If there are critical differences, it 
    will re-generate the parameters and update the pickle file.  This process 
    drastically improves the load time when running the same setup many times.
    To force the program to regenerate the parameters, provide the '-rp' flag to
    remove the pickle file.  This should occur when re-coarse-graining after a
    minimization.  This has not yet been automated.    
    '''
    pkl_file = ARGS.pdb[:-3] + 'pkl'

    if os.path.isfile(pkl_file) and ARGS.rm_pkl:
        print '>>> removing cg paramaters for %s: %s' % (ARGS.pdb, pkl_file)
        os.remove(pkl_file)

    do_cg_pro = do_cg_dna = do_dna_flex = False

    # if pickle_file exists:
    if os.path.isfile(pkl_file):
        print 'loading cg parameters for %s from: %s' % (ARGS.pdb, pkl_file)
        # load pickle_file
        pkl_in = open(pkl_file, 'rb')

        # critical
        cg_dna          = pickle.load(pkl_in)
        aa_dna          = pickle.load(pkl_in)
        cg_pro          = pickle.load(pkl_in)
        aa_pro          = pickle.load(pkl_in)
        vecXYZ          = pickle.load(pkl_in)

        trialbeads      = pickle.load(pkl_in)
        beadgroups      = pickle.load(pkl_in)

        move_masks      = pickle.load(pkl_in)
        all_beads       = pickle.load(pkl_in)
        dna_bead_masks  = pickle.load(pkl_in)
        aa_pgroup_masks = pickle.load(pkl_in)
        cg_pgroup_masks = pickle.load(pkl_in)
        all_proteins    = pickle.load(pkl_in)
        aa_all          = pickle.load(pkl_in)
        aa_pro_mask     = pickle.load(pkl_in)
        aa_dna_mask     = pickle.load(pkl_in)
        
        # input parametrs used to generate these cg-parameters
        ARGS_old        = pickle.load(pkl_in)        
        flex_resids_old = pickle.load(pkl_in)
        pro_groups_old  = pickle.load(pkl_in)
        
        pkl_in.close()
        
        # check if input parameters have changes since last using this pdb
        if pro_groups != pro_groups_old:
            do_cg_pro = True
            print '>>>Previous run parameter, %s, conflict with present run' % 'pro_groups'
            print '>>>Re-coarse-graining the Protein'               
            
        if ARGS.bp_per_bead != ARGS_old.bp_per_bead:
            do_cg_dna = True
            print '>>>Previous run parameter, %s, conflict with present run' % 'bp_per_bead'
            print '>>>Re-coarse-graining the DNA'       

        if flex_resids != flex_resids_old:
            do_dna_flex = True
            print '>>>Previous run parameter, %s, conflict with present run' % 'flex_resids'
            print '>>>Re-identifying flexible residues'       
            
    else:
        # load in the all atom pdb
        aa_all = sasmol.SasMol(0)
        aa_all.read_pdb(ARGS.pdb)  

        do_cg_dna = do_cg_pro = True
        

    if do_cg_dna:
        #~~~ DNA ONLY SECTION ~~~#
        (aa_dna, aa_dna_mask, cg_dna, all_beads, dna_bead_masks, 
         vecXYZ) = make_cg_dna(dna_segnames, dna_resids, ARGS.bp_per_bead, 
                               aa_all)
        do_dna_flex = True
        
    if do_dna_flex:
        #~~~ Determine flexible DNA beads from user input ~~~#
        beadgroups, trialbeads = is_bead_flexible(flex_resids, cg_dna.natoms(), 
                                    dna_resids[0], ARGS.bp_per_bead)

    if do_cg_pro:
        #~~~ Protein ONLY section ~~~#
        (aa_pro, aa_pro_mask, cg_pro, cg_pgroup_masks, aa_pgroup_masks, 
         all_proteins, move_masks) = make_cg_pro(aa_all, pro_groups)
        
    if do_cg_dna or do_cg_pro or do_dna_flex:
        print 'cg %s using updated parameters. Result saved to: %s' % (ARGS.pdb,
                                                                       pkl_file)

        # create the pickle_file for future use
        pkl_out = open(pkl_file, 'wb')

        pickle.dump(cg_dna, pkl_out, -1)
        pickle.dump(aa_dna, pkl_out, -1)
        pickle.dump(cg_pro, pkl_out, -1)
        pickle.dump(aa_pro, pkl_out, -1)
        pickle.dump(vecXYZ, pkl_out, -1)
        pickle.dump(trialbeads, pkl_out, -1)
        pickle.dump(beadgroups, pkl_out, -1) 
        pickle.dump(move_masks, pkl_out, -1) 
        pickle.dump(all_beads, pkl_out, -1) 
        pickle.dump(dna_bead_masks, pkl_out, -1)
        pickle.dump(aa_pgroup_masks, pkl_out, -1) 
        pickle.dump(cg_pgroup_masks, pkl_out, -1) 
        pickle.dump(all_proteins, pkl_out, -1)  
        pickle.dump(aa_all, pkl_out, -1)
        pickle.dump(aa_pro_mask, pkl_out, -1)    
        pickle.dump(aa_dna_mask, pkl_out, -1)    

        # input parameters used to get these cg-parameters
        pickle.dump(ARGS, pkl_out, -1)
        pickle.dump(flex_resids, pkl_out, -1)
        pickle.dump(pro_groups, pkl_out, -1)        
        
        pkl_out.close()
    
    # this is important for re-aligning the proteins after moving them
    cg_pro_orig = sasmol.SasMol(0)
    error = cg_pro.copy_molecule_using_mask(cg_pro_orig,
                                        np.ones(len(cg_pro.coor()[0])), 0)
    
    return (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, trialbeads, beadgroups, 
            move_masks, all_beads, dna_bead_masks, aa_pgroup_masks, 
            cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask)






if __name__ == "__main__":

    import argparse

    ''' wanted to implement this for error handling but never did'''
    # import logging
    # LOGGER = logging.getLogger(__name__) #add module name manually
    # import sys
    # if '-v' in sys.argv:
    #     logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
    # else:
    #     logging.basicConfig()

    ARGS = parse()  # this makes ARGS global

    main()
    print '\nFinished %d successful DNA moves! \n\m/ >.< \m/' % ARGS.nsteps
