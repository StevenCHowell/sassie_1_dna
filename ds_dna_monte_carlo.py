#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: generate modified DNA or DNA-protein structures
# Created: 1 Decemeber 2013

# $Id$
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sassie.sasmol.sasmol as sasmol
import sassie.simulate.monte_carlo.ds_dna_monte_carlo.collision as collision
import random, warnings, time, os, numpy
import sassie.simulate.monte_carlo.ds_dna_monte_carlo.special.input_filter as input_filter
from sassie.tools.merge_utilities import merge_dcd_files
import multiprocessing

try: 
    import cPickle as pickle
except:
    import pickle as pickle
    
# added for openmm
import sassie.sasconfig as sasconfig
# import sys,locale,subprocess
# from sys import stdout, exit, stderr
try:
    import sassie.simulate.openmm.openmm as omm
    from simtk.openmm.app import *
    from simtk.openmm import *
    import  simtk.unit as unit
except:
    print 'failed to load packages for running openmm'
    print 'BE AWARE: will continue to run with non-optimal setup'

'''
        DS_DNA_MONTE_CARLO is the module that performs Monte Carlo moves on DNA
        structures from a dcd/pdb.
'''

def print_failure(message,txtOutput):

        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
        txtOutput.put(message)

        return
    
def unpack_variables(variables):
    
    # standard user input
    try:
        ofile  = variables['ofile'][0]
    except:
        ofile  = variables['outfile'][0]
    try:
        infile = variables['infile'][0]
    except:
        infile = variables['pdbfile'][0]
    refpdb     = variables['refpdb'][0]
    path       = variables['path'][0]
    trials     = variables['trials'][0]
    goback     = variables['goback'][0]
    runname    = variables['runname'][0]
    psffile    = variables['psffile'][0]
    
    # molecule specific input    
    theta_max    = variables['theta_max'][0]
    theta_z_max  = variables['theta_z_max'][0]
    bp_per_bead  = variables['bp_per_bead'][0]
    dna_segnames = variables['dna_segnames'][0]
    dna_resids   = variables['dna_resids'][0]
    rigid_groups = variables['rigid_groups'][0]
    flex_resids  = variables['flex_resids'][0]
    
    # specialized/advanced input
    debug         = variables['debug'][0]
    write_flex    = variables['write_flex'][0]
    keep_cg_files = variables['keep_cg_files'][0]
    keep_unique   = variables['keep_unique'][0]
    n_dcd_write   = variables['n_dcd_write'][0]
    softrotation  = numpy.abs(variables['softrotation'][0])
    rm_pkl        = variables['rm_pkl'][0]
    openmm_min    = variables['openmm_min'][0]
    seed          = variables['seed'][0]
    temperature   = variables['seed'][0]

    return (ofile, infile, refpdb, path, trials, goback, runname, psffile, theta_max, 
            theta_z_max, bp_per_bead, dna_segnames, dna_resids, rigid_groups, 
            flex_resids, debug, write_flex, keep_cg_files, keep_unique,
            n_dcd_write, softrotation, rm_pkl, openmm_min, seed, temperature)

def make_cg_dna(dna_segnames, dna_resids, bp_per_bead, aa_all, txtOutput,
                dna_bead_masks=[],aa_dna_mask=[]):
    frame=0
    dna1 = dna_segnames[0]
    dna2 = dna_segnames[1]
    resid1 = dna_resids[0]
    resid2 = dna_resids[1]
    
    # check the input
    assert numpy.abs(resid1[1]-resid1[0]) == numpy.abs(resid2[1]-resid2[0]), (
        "number of paired bases in DNA strands are not the same")

    bps = numpy.abs(resid1[1]-resid1[0])+1
    nbeads = int(numpy.round(bps/float(bp_per_bead)))
    #print 'nbeads =', nbeads

    # load the dna into its own sasmol object
    aa_dna = sasmol.SasMol(0)
    if 0 == len(aa_dna_mask):
        dna_filter = '((moltype[i] == "dna" or moltype[i] == "rna"))' # RNA problem
        error, aa_dna_mask = aa_all.get_subset_mask(dna_filter)
    error = aa_all.copy_molecule_using_mask(aa_dna, aa_dna_mask, frame) # select dna
    natomsD = aa_dna.natoms()
    print 'dna atoms =', natomsD

    # populate the cg_dna sasmol object with semi-arbitrary atom properties
    cg_dna = sasmol.SasMol(0)
    basis_filter = ("((name[i] == 'N1') and (resid[i] >= "+str(resid1[0])+
                    " and resid[i] < "+str(resid1[0] + nbeads)+") and "
                    "(segname[i] == '"+dna1+"'))")
    error, mask = aa_dna.get_subset_mask(basis_filter)   #; print mask
    s = numpy.sum(mask)
    if s != nbeads:
        message = ("\n>>> ERROR!!! fail to create correct number of dummy"
                         " beads: expected %d but made %d\n") %(nbeads, s)
        print_failure(message, txtOutput)
    assert s == nbeads, ("\n>>> ERROR!!! fail to create correct number of dummy"
                         " beads: expected %d but made %d\n") %(nbeads, s)

    # initialize bead coordinates and other variables
    error = aa_dna.copy_molecule_using_mask(cg_dna, mask, frame)

    cg_coor = numpy.zeros((1, nbeads, 3))
    r1a = resid1[0]
    r2a = resid2[0]
    all_beads = []       # list of the all atom sasmol object for each bead
    # list of the all atom masks for each bead
    if dna_bead_masks:
        do_dna_bead_masks = False
    else:
        do_dna_bead_masks = True

    print "coarse-graining the DNA, this may take awhile..."
    tic = time.time()

    bp_per_bead = numpy.ones(nbeads,dtype=int) * bp_per_bead
    
    for j in xrange(nbeads):
        bead = sasmol.SasMol(0)
        if do_dna_bead_masks:
            if bp_per_bead[j] > 1:
                if j+1 == nbeads:
                    # accomodate for a non-divisible number of residues
                    bp_per_bead[j] = bps - bp_per_bead[j] * j 
    
                # Get the atoms from DNA strand 1
                if resid1[0] < resid1[1]:
                    r1b = r1a + bp_per_bead[j]  # r1b > r1a
                    basis_filter1 = ("((resid[i] > "+str(r1a-1)+" and "
                                     "resid[i] < "+str(r1b)+") and "
                                     "(segname[i]=='"+dna1+"')) or ")
                else:
                    r1b = r1a - bp_per_bead[j]  # r1b < r1a
                    basis_filter1 = ("((resid[i] > "+str(r1b)+" and "
                                     "resid[i] < "+str(r1a+1)+") and "
                                     "(segname[i]=='"+dna1+"')) or ")
        
                # Get the atoms from DNA strand 2
                if resid2[0] < resid2[1]:
                    r2b = r2a + bp_per_bead[j]  # r2b > r2a
                    basis_filter2 = ("((resid[i] > "+str(r2a-1)+" and "
                                     "resid[i] < "+str(r2b)+") and "
                                     "(segname[i]=='"+dna2+"'))")
                else:
                    r2b = r2a - bp_per_bead[j]   # r2b < r2a
                    basis_filter2 = ("((resid[i] > "+str(r2b)+" and "
                                     "resid[i] < "+str(r2a+1)+") and "
                                     "(segname[i]=='"+dna2+"'))")
            else:
                # Get the atoms from DNA strand 1
                if resid1[0] < resid1[1]:
                    r1b = r1a + bp_per_bead[j]  # r1b > r1a
                else:
                    r1b = r1a - bp_per_bead[j]  # r1b < r1a
                
                basis_filter1 = ("((resid[i] == "+str(r1a)+") and "
                                 "(segname[i]=='"+dna1+"')) or ")
        
                # Get the atoms from DNA strand 2
                if resid2[0] < resid2[1]:
                    r2b = r2a + bp_per_bead[j]  # r2b > r2a
                else:
                    r2b = r2a - bp_per_bead[j]   # r2b < r2a
    
                basis_filter2 = ("((resid[i]== "+str(r2a)+") and "
                                 "(segname[i]=='"+dna2+"'))")
                
            basis_filter = basis_filter1+basis_filter2
    
            # create a mask to select the atoms for the bead
            error, mask = aa_dna.get_subset_mask(basis_filter)
    
            # store the mask for the reverse coarse-graining
            dna_bead_masks.append(mask)

            # setup for next iteration
            r1a = r1b
            r2a = r2b
        else:
            mask = dna_bead_masks[j]
            
        error = aa_dna.copy_molecule_using_mask(bead, mask, frame)

        # not the best choice, COMs of bps in linear DNA isn't perfectly linear
        
        com = bead.calccom(0)
        cg_coor[0, j, :] = com  # a list of the com coordinate for each bead

        # new method
        # basis_filter = ("backbone")
        # error, backbone_mask = bead.get_subset_mask(basis_filter)
        # test_mol = sasmol.SasMol(0)
        # error = bead.copy_molecule_using_mask(test_mol, backbone_mask, 0)
        # com = test_mol.calccom(0)
        # cg_coor[0, j, :] = com

        # calculate atomic coodinates, using the bead com as the origin
        bead.center(0)
        all_beads.append(bead)
        
    cg_dna.setCoor(cg_coor)  #; print cg_dna._coor[frame, :, 0]

    vecXYZ = numpy.zeros((3, nbeads, 3))
    vecXYZ[0] = [1, 0, 0]
    vecXYZ[1] = [0, 1, 0]
    vecXYZ[2] = [0, 0, 1]
    toc = time.time() - tic
    print 'DNA coarse-graining took %0.3f seconds' % toc
    return (aa_dna, aa_dna_mask, cg_dna, all_beads, dna_bead_masks, vecXYZ,
            bp_per_bead)

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
            cg_pgroup_mask = numpy.zeros(sum(aa2cg_pro_mask), dtype='int32')
            aa_pgroup_mask = numpy.zeros(len(aa2cg_pro_mask), dtype='int32')

        cg_pgroup_masks.append(cg_pgroup_mask)
        aa_pgroup_masks.append(aa_pgroup_mask)
        all_proteins.append(aa_pro_group)

    # combined the protein groups into move groups
    move_masks = []
    for i in xrange(len(cg_pgroup_masks)):
        move_masks.append(numpy.copy(cg_pgroup_masks[i]))
    for i in xrange(len(move_masks)):
        for j in xrange(i+1, len(move_masks)):
            move_masks[i] += move_masks[j]
            
    toc = time.time() - tic
    print 'Protein coarse-graining took %0.3f seconds' % toc
    
    return (aa_pro, aa_pro_mask, cg_pro, cg_pgroup_masks, aa_pgroup_masks, 
            all_proteins, move_masks)

def make_rigid_groups(aa_all, rigid_groups):
    # group the rigid bodies
    
    print "grouping the rigid bodies..."
    tic = time.time()
    if 0 < len(rigid_groups):
        # create a sasmol object for each rigid move group (rigid_groups)
        rigid_group_masks = []  # masks to separate the cg_protein into NCP groups
        rigid_group_mols = []
        (pre, post, btwn) = ("((segname[i]=='", "'))", " or ")
        for group in rigid_groups:
            rigid_group_mol = sasmol.SasMol(0)
            if len(group) > 0:
                basis_filter = ""
                for seg in group:
                    basis_filter += pre + seg + post + btwn
                basis_filter = basis_filter[0:-4]  # remove the last btwn
                error, rigid_group_mask = aa_all.get_subset_mask(basis_filter)
                error = aa_all.copy_molecule_using_mask(rigid_group_mol, 
                                                        rigid_group_mask, 0)
            else:
                rigid_group_mask = numpy.zeros(aa_all.natoms(), dtype='int32')
    
            rigid_group_masks.append(rigid_group_mask)
            rigid_group_mols.append(rigid_group_mol)
    
        # combined the protein groups into move groups
        rigid_move_masks = []
        for i in xrange(len(rigid_group_masks)):
            rigid_move_masks.append(numpy.copy(rigid_group_masks[i]))
        for i in xrange(len(rigid_move_masks)):
            for j in xrange(i+1, len(rigid_move_masks)):
                rigid_move_masks[i] += rigid_move_masks[j]
                
        rigid_mask = numpy.copy(rigid_move_masks[0])
    else:
        rigid_mask = numpy.zeros(aa_all.natoms())
        rigid_group_masks = []
        rigid_group_mols = []
        rigid_move_masks = []
        
    rigid_mol = sasmol.SasMol(0)
    error = aa_all.copy_molecule_using_mask(rigid_mol,rigid_mask, 0)
    toc = time.time() - tic
    print 'Grouping the rigid bodies took %0.3f seconds' % toc
    
    return (rigid_mol, rigid_mask, rigid_group_masks, rigid_group_mols, 
            rigid_move_masks)

def is_bead_flexible(flex_resids, nbeads, resid1, bp_per_bead, debug=False):
    ## Flexible part
    # create an array (groupid) labeling the group for each flexible cg-DNA bead
    flexResids = [item for sublist in flex_resids for item in sublist]
    groupid = numpy.zeros(len(flexResids))
    link = numpy.zeros((nbeads-1, 2), dtype=int)  # links between DNA beads
    # if a bead is designated flexible both it's links are flexible

    n_start = 0
    for i in xrange(len(flex_resids)):
        n_end = len(flex_resids[i]) + n_start
        groupid[n_start:n_end] = numpy.ones(len(flex_resids[i])) * i
        n_start = n_end
    
    flexResids = numpy.concatenate(([flexResids], [groupid]),0).T    

    r1a = resid1[0]
    for j in xrange(nbeads):
        if resid1[0] < resid1[1]:
            r1b = r1a + bp_per_bead[j]  # r1b > r1a
        else:
            r1b = r1a - bp_per_bead[j]  # r1b < r1a            
        ## check residues from DNA strand 1 to see if this is a flexible bead
        for i in xrange(r1a, r1b):
            if i in flexResids[:, 0]:
                (r, c) = numpy.where(flexResids==i)
                # this stores the assigned group in the temp var: beadGroup
                beadGroup = flexResids[r[0], 1]
                if j==0:                         # --> first bead
                    link[j] = [1, beadGroup]     # first link flexible
                elif j==nbeads-1:                # --> last bead
                    link[-1] = [1, beadGroup]    # last link flexible
                else:                            # --> all other beads
                    link[j-1:j+1] = [1, beadGroup]
                break  # bead defined as flexible, stop looping over residues        

        # setup for next iteration
        r1a = r1b

    # position the bead at its calculated com
    beadgroups = numpy.zeros(len(link[:, 1])+1, dtype=int)
    beadgroups[1:] = link[:, 1]
    # need a trial bead for every flexible link (0th bead will never move)
    trialbeads = numpy.zeros(numpy.sum(link[:, 0]), dtype=numpy.uint)
    j = 0
    for i_link in xrange(nbeads-1):
        if link[i_link, 0]:
            trialbeads[j] = i_link+1
            j+=1

    if debug:
        print 'beadgroups =', beadgroups[trialbeads]

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
    cg_natoms = numpy.copy(cg_dna.natoms())
    coor = numpy.copy(cg_dna.coor()[0])

    # split vecXYZ into three matrices
    vecX = numpy.copy(vecXYZ[0])
    vecY = numpy.copy(vecXYZ[1])
    vecZ = numpy.copy(vecXYZ[2])

    for i in xrange(cg_natoms):
        # R would align the rotated coordinates back to where they originated
        R = align2xyz(vecX[i, :], vecY[i, :], vecZ[i, :])

        # M will rotate cooridates from their original reference frame to the
        # rotated refrenece frame of the bead
        # then translate the com to the bead's com
        M = R.transpose()
        M[3, :3] = coor[i, :] # insert the com translation into the M matrix

        r, c = allBeads[i].coor()[0].shape
        beadCoor = numpy.ones((r, 4))           # initialize the beadCoor matrix
        origCoor = numpy.ones((1, r, 3))
        origCoor = allBeads[i].coor()
        beadCoor[:, :3] = numpy.copy(origCoor[0])
        beadCoor = numpy.dot(beadCoor, M)

        tmpCoor = numpy.zeros((1, r, 3))        # initialize the tmpCoor
        tmpCoor[0] = beadCoor[:, :3]
        allBeads[i].setCoor(tmpCoor)

        #recombine all the beads back into one pdb
        e = aa_dna.set_coor_using_mask(allBeads[i], 0, masks[i])
        error.append(e)

        #reset the bead coordinates for the next iteration
        allBeads[i].setCoor(origCoor)

    return error

def align2xyz(vecX, vecY, vecZ):

    zero = 10e-5
    
    tmp_coor = numpy.zeros((2, 4))
    tmp_coor[:, 3] = 1
    tmp_coor[1, 0:3] = vecZ
    A1 = align2z(tmp_coor)

    newX = numpy.dot(vecX, A1[0:3, 0:3]) #; print 'newX =', newX
    #newY = numpy.dot(vecY, A1[0:3, 0:3]) ;# print 'newY =', newY
    if newX[2] >= zero:
        message = "ERROR!!! z-component of newX is not zero and it should be"
        print_failure(message, txtOutput)
        assert newX[2] < zero, ("ERROR!!! z-component of newX is not zero and it "
                             "should be")

    thetaZ_x = -numpy.arctan2(newX[1], newX[0])
    #thetaZ_y = -numpy.arctan2(newY[1], -newY[0])

    A2 = rotate4x4('z', thetaZ_x)  #; print 'A2 =', A2

    A = numpy.dot(A1, A2)

    newY = numpy.dot(vecY, A[0:3, 0:3]) #; print 'finalY =', newY
    if newY[0]+newY[2] > 1 + zero:
        message = ("ERROR!!! newY is not aligned to the "
                                       "y-axis and it should be")
        print_failure(message, txtOutput)
        assert newY[0]+newY[2] < 1 + zero, ("ERROR!!! newY is not aligned to the "
                                       "y-axis and it should be")

    return A

def rotate4x4(axis, theta):
    R = numpy.eye(4)
    ct = numpy.cos(theta)
    st = numpy.sin(theta)
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

    # print 'R =', R

    return R

def move2origin(coor4):
    T    = numpy.eye(4, dtype=numpy.float)
    Ti   = numpy.copy(T)
    #s print 'coor passed to T, Ti:\n', coor4
    #s print 'size =', coor4.shape
    if not coor4.shape < (1, 4):
        T[3, 0:3] = -coor4[0, 0:3]
        Ti[3, 0:3] = coor4[0, 0:3]

    #s print 'T:\n', T
    #s print 'test:\n', numpy.dot(coor4, T)
    #s print 'Ti:\n', Ti
    #s print 'testI:\n', numpy.dot(numpy.dot(coor4, T), Ti)
    return (T, Ti)

def align2z(coor4):
    '''
    align the axis connecting the first 2 coordinates of a 1x4 array of
    coodinate vectors to the z-axis
    '''
    A = numpy.eye(4, dtype=numpy.float)
    #s Axz  = numpy.eye(4, dtype=numpy.float)
    #s Az   = numpy.eye(4, dtype=numpy.float)
    assert all(coor4[0] == [0., 0., 0., 1., ]), ("coordinates passed to align2z"
                    "were not translated to the origin")

    if coor4.shape > (1, 4):
        (u, v, w) = coor4[1, 0:3]
        small = 1E-14 # small limit to make sure not dividing by zero

        d1 = numpy.sqrt(u**2+v**2)
        d2 = numpy.sqrt(u**2+v**2+w**2)
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
        #s A = numpy.dot(Axz, Az) #;   print 'A3=\n', A
    return A

def beadRotate(coor3, vecXYZ, thetas, rigid_coor3, n_soft=1):
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
    coor4 = numpy.ones((natoms, 4), numpy.float)
    X = numpy.copy(coor4)
    Y = numpy.copy(coor4)
    Z = numpy.copy(coor4)
    coor4[:, 0:3] = coor3     #; print 'coor4 =', coor4
    X[:, 0:3] = numpy.copy(vecXYZ[0])  #;print 'X = \n', X
    Y[:, 0:3] = numpy.copy(vecXYZ[1])
    Z[:, 0:3] = numpy.copy(vecXYZ[2])

    rigid_coor4 = numpy.ones((len(rigid_coor3), 4), numpy.float)
    rigid_coor4[:, 0:3] = rigid_coor3

    # create the translation-rotation matrix
    # This is intended to be multiplied from the right (unlike standard matrix
    # multiplication) so as not to require transposing the coordinate vectors.
    # print '(tx, ty, tz) in degrees:', thetas
    tx = thetas[0]*(numpy.pi/180.0)
    ty = thetas[1]*(numpy.pi/180.0)
    tz = thetas[2]*(numpy.pi/180.0)
    # print '(tx, ty, tz) in radians:', (tx, ty, tz)
    cx = numpy.cos(tx)
    sx = numpy.sin(tx)
    cy = numpy.cos(ty)
    sy = numpy.sin(ty)
    cz = numpy.cos(tz)
    sz = numpy.sin(tz)

    # initialize the rotation
    # consolidated method of defining the rotation matrices
    R = numpy.eye(4, dtype=numpy.float)
    (R[0][0], R[0][1], R[0][2]) = ( cy*cz, cy*sz, -sy )
    (R[1][0], R[1][1], R[1][2]) = ( sx*sy*cz-cx*sz, sx*sy*sz+cx*cz, sx*cy )
    (R[2][0], R[2][1], R[2][2]) = ( sx*sz+cx*sy*cz, cx*sy*sz-sx*cz, cx*cy )

    # verbose method of defining the rotation matrices
    #s Rx = numpy.eye(4, dtype=numpy.float)
    #s Ry = numpy.eye(4, dtype=numpy.float)
    #s Rz = numpy.eye(4, dtype=numpy.float)
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
    #s R = numpy.dot(numpy.dot(Rx, Ry), Rz) #; print "Rxyz = \n", Rxyz

    (T0, Ti0) = move2origin(coor4)  # Ti0 != T0, negative off diag elements
    coor4 = numpy.dot(coor4, T0)  #; print 'moved to origin coor:\n', coor4
    A = align2z(coor4)
    AR = numpy.dot(A, R)

    # first of n_soft rotations
    coor4 = numpy.dot(coor4, AR)

    # the coarse grained beads local coordinates should not be translated,
    # only rotated
    X = numpy.dot(X, AR)
    Y = numpy.dot(Y, AR)
    Z = numpy.dot(Z, AR)

    rigid_coor4 = numpy.dot(numpy.dot(rigid_coor4, T0), AR)

    # remaining of n_soft rotations
    for i in xrange(1, n_soft):
        (T, Ti) = move2origin(coor4[i:])
        coor4[i:] = numpy.dot(coor4[i:], T)  # move to origin
        # print 'moved2origin coor:\n', coor4
        rigid_coor4 = numpy.dot(rigid_coor4, T)

        coor4[i:] = numpy.dot(coor4[i:], R)
        # print 'step %d' %i, 'rotated coor:\n', coor4
        rigid_coor4 = numpy.dot(rigid_coor4, R)

        coor4[i:] = numpy.dot(coor4[i:], Ti) # return to original position
        # print 'returned from origin coor:\n', coor4
        rigid_coor4 = numpy.dot(rigid_coor4, Ti)

        # coarse grained beads local coordinates should not be translated,
        # only rotated
        X[i:] = numpy.dot(X[i:], R)
        Y[i:] = numpy.dot(Y[i:], R)
        Z[i:] = numpy.dot(Z[i:], R)

    coor4 = numpy.dot(coor4, A.transpose())
    coor4 = numpy.dot(coor4, Ti0)
    X = numpy.dot(X, A.transpose())
    Y = numpy.dot(Y, A.transpose())
    Z = numpy.dot(Z, A.transpose())

    rigid_coor4 = numpy.dot(rigid_coor4, A.transpose())
    rigid_coor4 = numpy.dot(rigid_coor4, Ti0)

    vecXYZ[0, 1:] = X[1:, 0:3]
    vecXYZ[1, 1:] = Y[1:, 0:3]
    vecXYZ[2, 1:] = Z[1:, 0:3]


    # this returns the modified positions and orientations for all but the
    # first bead or reference bead
    return (coor4[1:, 0:3], vecXYZ[:, 1:], rigid_coor4[:, 0:3])

def checkU(coor):
    '''
    method to generate and check the u-vectors pointing from bead to bead
    (used to make sure they were are all the same length but this did not work
    with real coarse-grained DNA)
    returns the u-vectors and their average magnitude (l)
    '''
    u = coor[1:, :]-coor[:-1, :]                  # u-vectors between beads
    lu = numpy.sqrt(u[:, 0]**2+u[:, 1]**2+u[:, 2]**2) # magnitues of u-vectors
    nu = lu.size
    u = u/lu.reshape(nu,1)   # this make the u-vectors into unit vectors

    # Commented this out because it did not work for real DNA because base pairs
    # did not always start evenly spaced apart
    #s erLim = 1e-3
    #s test = numpy.abs(lu-l) > erLim  # check if any of the l-lengths are different
    #s if test.any():
    #s         print 'ERROR: the beads are not uniformly spaced'
    #s         print 'u = \n', u
    #s         print test
    #s         # print 'difference =', u[:, 2] - l
    #s         print 'distance =', lu

    return (u, lu)

def energyBend(lp, u, l):
    '''
    this function finds the E/kT for N beads connected by N-1 rods
    lp is the persistence length
    u contains the N-1 unit vectors representing each rod
    l contains the distances between beads
    ### this has the correted bending energy ###
    '''

    uk0 = u[:-1, :]          #; print 'u_k= \n', uk0
    uk1 = u[1:, :]             #; print 'u_k= \n', uk1
    res = numpy.sum((1-numpy.inner(uk0, uk1).diagonal())/l[:-1])  # calculate the sum of their products

    return res*lp

def f_energy_wca(w, coor, wca0, trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN
    '''

    wca1 = numpy.copy(wca0)
    # calculate the distance between moved DNA beads and all other DNA beads
    # increment trial bead by one because it will index from 1 instead of 0
    wca1 = collision.wca_d(coor, trial_bead+1, w, wca1)

    res = 4.*numpy.sum(wca1)
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

    mask_indices = numpy.nonzero(mask*numpy.arange(1, len(mask)+1))[0]

    return mask_indices

def dna_mc(trials, i_loop, theta_max, theta_z_max, debug, goback, n_dcd_write, 
           keep_unique, keep_cg_files, softrotation, write_flex, runname, ofile,
           cg_dna, aa_dna, rigid_mol, vecXYZ, lp, trialbeads, beadgroups, 
           rigid_move_masks, all_beads, dna_bead_masks, rigid_group_masks, 
           rigid_group_mols, aa_all, rigid_mask, aa_dna_mask, seed=0, 
           dna_type='b'):
    '''
    this function perform nsteps Monte-Carlo moves on the cg_dna
    '''

    if softrotation > 1 and debug:
        print 'softening all rotaitons over %d beads' % softrotation

    # timestr = time.strftime("_%y%m%d_%H%M%S") # suffix for output files

    if debug and False:
        # never got this to work right
        print 'sending coordinates to vmd port 2222'
        aa_all.send_coordinates_to_vmd(2222,0)
    
    dna_path = runname + '/dna_mc/'
    direxist = os.path.exists(dna_path)
    if(direxist==0):
        os.system('mkdir -p ' + dna_path)
        
    print 'runname =', runname

    if write_flex:
        write_flex_resids(infile, all_beads, trialbeads, dna_path)

    aa_ofile = dna_path + ofile[:-4] + '_%03d.dcd' % i_loop
    aa_all_dcd_out = aa_all.open_dcd_write(aa_ofile)
    
    # create the coarse-grained DNA and protein dcd and pdb files
    cg_dna_ofile = dna_path + 'cg_dna' + '_%03d.dcd' % i_loop
    rigid_ofile = dna_path + 'rigid' + '_%03d.dcd'% i_loop
    cg_dna.write_pdb(cg_dna_ofile[:-4] + '.pdb', 0, 'w')
    rigid_mol.write_pdb(rigid_ofile[:-4] + '.pdb', 0, 'w')    
    cg_dna_dcd_out = cg_dna.open_dcd_write(cg_dna_ofile)
    rigid_dcd_out = rigid_mol.open_dcd_write(rigid_ofile)    
    
    # will write these out to dcd files to store the coordinates along the way
    vecX_mol = sasmol.SasMol(0)
    vecY_mol = sasmol.SasMol(0)
    vecZ_mol = sasmol.SasMol(0)    
    error, mask = cg_dna.get_subset_mask('(all)')
    error = cg_dna.copy_molecule_using_mask(vecX_mol,mask,0)
    error = cg_dna.copy_molecule_using_mask(vecY_mol,mask,0)
    error = cg_dna.copy_molecule_using_mask(vecZ_mol,mask,0)
    vecX_mol.setCoor(numpy.array([vecXYZ[0]])) # the numpy.array recast these so they 
    vecY_mol.setCoor(numpy.array([vecXYZ[1]])) # do not update with vecXYZ 
    vecZ_mol.setCoor(numpy.array([vecXYZ[2]]))
    vecX_dcd_name = dna_path + 'vecX' + '%03d.dcd' % i_loop
    vecY_dcd_name = dna_path + 'vecY' + '%03d.dcd' % i_loop
    vecZ_dcd_name = dna_path + 'vecZ' + '%03d.dcd' % i_loop
    vecX_dcd_out = vecX_mol.open_dcd_write(vecX_dcd_name)
    vecY_dcd_out = vecY_mol.open_dcd_write(vecY_dcd_name)    
    vecZ_dcd_out = vecZ_mol.open_dcd_write(vecZ_dcd_name)        
    vecX_mol.write_dcd_step(vecX_dcd_out, 0, 1)
    vecY_mol.write_dcd_step(vecY_dcd_out, 0, 1)   
    vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, 1)   
    
    # initialize variables for each run
    steps_from_0 = numpy.zeros(trials, dtype='int64')
    xyz = numpy.copy(vecXYZ)
    d_coor = numpy.copy(cg_dna.coor()[0]) # unique memory for each
    r_coor = numpy.copy(rigid_mol.coor()[0]) # unique memory for each

    (u, l) = checkU(d_coor) # vectors, u, and distances, l, between beads
    #s print "(u, l) =", (u, l) # debug info

    # effective width of the DNA chain
    # 4.6nm for B-Form (taken from: Odijk, T. Phys. Rev. E 2008, 77, 060901(R))
    dna_energy_width = {'a': 0, 'b': 46., 'z': 0} # yet to use a, and z type dna
    w = dna_energy_width[dna_type.lower()]
    if w > l.any():
        w = numpy.floor(l)
        if debug:
            print '>>> chain width, w, set to %d Angstroms (so w is < dist btwn beads)' %w

    dna_bead_radius = 4.5

    #!# this needs to be updated for generic overlap checking
    pro_bead_radius = 1.0 # 2A min seperation of CA atoms in database
    rigid_radius = 1.0
    
    pro_pro_test = pro_bead_radius + pro_bead_radius
    dna_pro_test = dna_bead_radius + pro_bead_radius
    rigid_rigid_test = rigid_radius + rigid_radius
    dna_rigid_test = dna_bead_radius + rigid_radius
    
    
    # calculate the energy of the starting positions
    wca0 = numpy.zeros((cg_dna.natoms(),cg_dna.natoms()))
    Ub0 = energyBend(lp, u, l)

    (Uwca0, wca0) = f_energy_wca(w, d_coor, wca0, 0)

    U_T0 = Ub0 + Uwca0
    # print '(Ub0, Uwca0, Ub0/U_T0, Uwca0/U_T0) =', (Ub0, Uwca0, Ub0/U_T0, 
    #                                                 Uwca0/U_T0)
    n_accept   = 0 # total times configuration was accepted
    n_reject   = 0 # total times configuration was rejected
    n_written  = 0 # total times dcd write has been called
    fail_tally = 0 # number of times failed for particular iteration
    n_from_reload = 0 # number of stps since last reload
    n_reload = [0]  # listt containing the i_goback values
    
    # this should not actually be >=, come back to this
    assert numpy.size(theta_max) - 1 >= numpy.max(beadgroups), (
        'each group needs its own theta_max: %d < %d'
        % (numpy.size(theta_max) - 1, numpy.max(beadgroups) ))

    # Main MC loop #
    while n_accept < trials:

        # Choose a bead to rotate
        trialbead = trialbeads[int((trialbeads.size)*numpy.random.random())]        

        # Determine rotation to perform
        trial_theta_max = theta_max[beadgroups[trialbead]]
        trial_theta_z_max = theta_z_max[beadgroups[trialbead]]
        thetaZ = 2 * trial_theta_z_max * numpy.random.random() - trial_theta_z_max
        thetaX = 2 * trial_theta_max   * numpy.random.random() - trial_theta_max
        thetaY = 2 * trial_theta_max   * numpy.random.random() - trial_theta_max
        thetaXYZ = [thetaX/softrotation,thetaY/softrotation,thetaZ/softrotation]

        if  len(rigid_move_masks) == 0 or beadgroups[trialbead] == len(rigid_move_masks):
            # Only DNA will be moving, create place-holder dummy coordinates
            r_coor_rot = numpy.zeros((0, 3))
        else:
            group = beadgroups[trialbead]
            r_mask_rot = rigid_move_masks[group]
            r_ind_rot = mask2ind(r_mask_rot)
            r_coor_rot = r_coor[r_ind_rot]

            if beadgroups[trialbead] == 0:
                # none of rigid componets will be fixed
                r_mask_fix = numpy.zeros(r_mask_rot.shape)
            else:
                r_mask_fix = rigid_move_masks[group-1] - rigid_move_masks[group]
            r_ind_fix = mask2ind(r_mask_fix)
            r_coor_fix = r_coor[r_ind_fix]

        # generate a newly rotated model
        (d_coor[trialbead:], xyz[:, trialbead:], r_coor_rot) = beadRotate(
            d_coor[trialbead-1:], xyz[:, trialbead-1:], thetaXYZ,
            r_coor_rot, softrotation) 

        # store the rotated protein coordinates
        if beadgroups[trialbead] < len(rigid_move_masks):
            r_coor[r_ind_rot] = r_coor_rot

        # calculate the change in energy (dU) and the boltzman factor (p)
        (u, l) = checkU(d_coor)
        Ub1 = energyBend(lp, u, l)

        # ~~~~ DNA interaction energy  ~~~~~~#
        (Uwca1, wca1) = f_energy_wca(w, d_coor, wca0, trialbead)

        U_T1 =  Ub1 + Uwca1
        dU = U_T1 - U_T0

        with warnings.catch_warnings():
            warnings.filterwarnings('error') # need this for np warnings
            try:
                probability = numpy.exp(-dU)
            # print '\n(Ub1, Uwca1) =', (Ub1, Uwca1) 
            # print '(Ub1/U_T1, Uwca1/U_T1) =', (Ub1/U_T1, Uwca1/U_T1)
            # print '(p, dU) =', (p, dU)
            except Warning:
                if dU > 99:
                    probability =  0
                    #s print 'energy was large, setting probability to 0'
                elif dU < 0:
                    probability =  1
                    #s print 'energy was negative, setting probability to 1'
                else:
                    print 'Warning: ~~> unclear OverflowError <~~ dU =', dU
                    print 'not sure where the error originated from'

        test = numpy.random.random()
        collision = 0

        if test >= probability:
            dna_pass = False
            # print 'step failed because of DNA energy'
        else:
            dna_pass = True

            # now check for collisions
            if len(r_coor_rot) > 0:   # only if proteins were rotated
                # ~~~~ Check for overlap, DNA-protein or protein-protein ~~~~~~#
                d_coor_fix = d_coor[trialbead:]
                d_coor_rot = d_coor[:trialbead]
                
                # check for protein-protein overlap
                if 1 == f_overlap2(r_coor_rot, r_coor_fix, rigid_rigid_test):
                    print 'Collision between 2 rigid components'
                    collision = 1
                    
                # check for DNA-protein overlap
                elif 1 == f_overlap2(r_coor_rot, d_coor_fix, dna_rigid_test):
                    print 'Rigid-DNA (rot-fix) collision'
                    collision = 1
    
                elif 1 == f_overlap2(r_coor_fix, d_coor_rot, dna_rigid_test):
                    print 'Rigid-DNA (fix-rot) collision'
                    collision = 1
    
        if dna_pass and collision == 0:
            n_from_reload += 1
            steps_from_0[n_accept] = n_from_reload + n_reload[-1]
            n_accept += 1                      # increment accept counter
            # cg_dna.setCoor(d_coor) # <-- DO NOT use setCoor, want uniuqe mem
            # cg_pro.setCoor(p_coor) # <-- DO NOT use setCoor, want uniuqe mem            
            cg_dna.setCoor(numpy.array([d_coor])) # update dna coordinates
            rigid_mol.setCoor(numpy.array([r_coor])) # update protein coordinates
            vecXYZ = numpy.copy(xyz)              # update dna orientations
            vecX_mol.setCoor(numpy.array([vecXYZ[0]])) # independent of vecXYZ[0]
            vecY_mol.setCoor(numpy.array([vecXYZ[1]])) # independent of vecXYZ[1]
            vecZ_mol.setCoor(numpy.array([vecXYZ[2]])) # independent of vecXYZ[2]
            
            wca0 = numpy.copy(wca1)               # update DNA WCA energy        
            U_T0 = U_T1                        # update total energy

            if debug:
                #print output regarding trial
                print "trial_bead(%3d) = %2d\t failed attempts = %2d" % (
                    n_accept, trialbead, fail_tally)
            else:
                print '.', ;
                
            fail_tally = 0                     # reset fail_tally
            
            # recover an all atom representation and save coordinates to a dcd
            # this requires re-inserting the aa-coordinates which takes added 
            # time so only do when designated
            if 0 == n_accept % n_dcd_write:
                # ~~recover aa-DNA~~
                error = recover_aaDNA_model(cg_dna, aa_dna, vecXYZ, all_beads,
                                            dna_bead_masks)
                # ~~recover aa-Protein~~
                # recover_aaPro_model(aa_pgroup_masks, rigid_group_masks, rigid_mol,
                                    # rigid_group_mols, aa_pro)
                                    
                # ~~Combine aa Complete Structure~~
                aa_all.set_coor_using_mask(rigid_mol, 0, rigid_mask)
                aa_all.set_coor_using_mask(aa_dna, 0, aa_dna_mask)
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written)
    
                # write out the accepted configuration for go-back use
                if goback > 0:
                    # default goback is -1 so this returns FALSE without user input
    
                    cg_dna.write_dcd_step(cg_dna_dcd_out, 0, n_written)
                    rigid_mol.write_dcd_step(rigid_dcd_out, 0, n_written)
                    # these are incremented by one because the 0th contains the 
                    # original coordinates 
                    vecX_mol.write_dcd_step(vecX_dcd_out, 0, n_written+1)
                    vecY_mol.write_dcd_step(vecY_dcd_out, 0, n_written+1)   
                    vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, n_written+1)                   

        else :
            if fail_tally == goback:  
                i_goback = rewind(seed, n_accept, cg_dna_ofile,
                            cg_dna, rigid_ofile, rigid_mol, vecX_dcd_name, 
                            vecX_mol, vecY_mol, vecY_dcd_name, vecZ_mol, 
                            vecZ_dcd_name, vecXYZ)

                d_coor = numpy.copy(cg_dna.coor()[0]) # reset the dna coordinates
                
                # reset the reference energy
                (u, l) = checkU(d_coor) 
                Ub0 = energyBend(lp, u, l)
                (Uwca0, wca0) = f_energy_wca(w, d_coor, wca0, 0)
                U_T0 =  Ub0 + Uwca0
                
                n_from_reload = 0
                n_reload.append(steps_from_0[i_goback-1])
                fail_tally = 0 # reset the fail counter
            else:
                fail_tally += 1                 # increment bead reject counter 
                n_reject += 1                   # increment total reject counter
                d_coor = numpy.copy(cg_dna.coor()[0]) # reset the dna coordinates
                
            r_coor = numpy.copy(rigid_mol.coor()[0]) # reset the protein coordinates
            xyz = numpy.copy(vecXYZ)              # reset the dna orientations

            # save previous coordinates again
            if not keep_unique:
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written )

    aa_all.close_dcd_write(aa_all_dcd_out)
    
    os.remove(vecX_dcd_name) 
    os.remove(vecY_dcd_name) 
    os.remove(vecZ_dcd_name) 

    if keep_cg_files:
        cg_dna.close_dcd_write(cg_dna_dcd_out)
        if len(rigid_move_masks) > 0:
            rigid_mol.close_dcd_write(rigid_dcd_out)
        else:
            os.remove(rigid_ofile[:-4] + '.pdb')
            os.remove(rigid_ofile)
    else:
        os.remove(cg_dna_ofile[:-4] + '.pdb')
        os.remove(cg_dna_ofile)
        os.remove(rigid_ofile[:-4] + '.pdb')
        os.remove(rigid_ofile)

    if goback > 0 and debug:
        goback_ofile = dna_path + 'n_from_0_%03d.txt' % i_loop
        numpy.savetxt(goback_ofile, steps_from_0, fmt='%d')

    if debug:
        print "accepted %d moves" % n_accept
        print "rejected %d moves" % n_reject

    print "accepted %0.2f percent of trials" % (100.0 * n_accept/(n_reject + n_accept))
    # print n_reload
    # print steps_from_0
    
    return aa_ofile, cg_dna_ofile, rigid_ofile

def rewind(seed, n_accept, cg_dna_ofile, cg_dna, cg_pro_ofile, cg_pro, 
           vecX_dcd_name, vecX_mol, vecY_mol, vecY_dcd_name, vecZ_mol,
           vecZ_dcd_name, vecXYZ):
    
    if seed > 0:
        numpy.random.seed(seed)
        i_goback = int(n_accept*numpy.random.random())
    else:
        i_goback = int(n_accept*numpy.random.random())
        
    # i_goback = 0 represents the original structure this is the 1st frame of 
    # the cg_dna_dcd and cg_pro_dcd so increment by 1
        
    vecX_mol.read_single_dcd_step(vecX_dcd_name, i_goback + 1)
    vecY_mol.read_single_dcd_step(vecY_dcd_name, i_goback + 1)   
    vecZ_mol.read_single_dcd_step(vecZ_dcd_name, i_goback + 1)                                     
    vecXYZ[0] = vecX_mol.coor()[0]
    vecXYZ[1] = vecY_mol.coor()[0]
    vecXYZ[2] = vecZ_mol.coor()[0]                

    if i_goback == 0:
        cg_dna.read_pdb(cg_dna_ofile[:-4] + '.pdb')
        cg_pro.read_pdb(cg_pro_ofile[:-4] + '.pdb')
        print '\n~~~ reloaded original coordinates ~~~\n'
    else:	
        cg_dna.read_single_dcd_step(cg_dna_ofile, i_goback + 1)
        cg_pro.read_single_dcd_step(cg_pro_ofile, i_goback + 1) 
        print '\n~~~ reloaded accepted coordinates #%d ~~~\n' % i_goback
        
    return i_goback

def write_flex_resids(infile, all_beads, flex_beads, dna_path=''):
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

    bases_per_bead = numpy.zeros((len(flex_beads)), dtype='int')
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

    flex_id_out = open(dna_path + infile[:-3]+'flex', 'w')
    flex_id_out.write('%s %s\n' % (dna1_label[0], dna2_label[0]))

    n_flex_bp = len(flex_beads)*bases_per_bead[0] / 2
    dna1_resids = numpy.reshape(numpy.array(dna1_bead_resids), n_flex_bp)
    dna2_resids_temp = numpy.reshape(numpy.array(dna2_bead_resids), n_flex_bp)
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
    flex_resids = numpy.genfromtxt(flex_file, dtype='int',delimiter=" ")[1:]
    return (segnames, flex_resids)

def main(variables):

    txtOutput=multiprocessing.JoinableQueue()

    (ofile, infile, refpdb, path, trials, goback, runname, psffile, theta_max, 
     theta_z_max, bp_per_bead, dna_segnames, dna_resids, rigid_groups, 
     flex_resids, debug, write_flex, keep_cg_files, keep_unique,
     n_dcd_write, softrotation, rm_pkl, openmm_min, seed, temperature
     ) = unpack_variables(variables)
    
    
    # set the DNA properties parameters:
    lp = 530.     # persistence length  (lp = 530A)
    # defined in dna_mc
    #w = 46        # effective measure of dsDNA chain width in A (w = 46A)
    #DOI 10.1021/ma201277e used w=46, lp=530

    if keep_unique:
        print 'only keeping unique structures'
    else:
        print 'keeping all accepted structures'
        
    # this is currently non-functional
    # if ARGS.Llp:
    #     print 'simulating long DNA'
    #     Llp = ARGS.Llp    # Llp is the length in units of lp: Llp = L/lp
    #     L = Llp*lp  #
    #     (cg_dna, vecXYZ) = makeLongDNA(Llp) # use this to make long cgDNA
    #     all_atom_pdb = '%d_llp' % Llp

    
    dna_resids =  []
    rigid_groups =  []

    if  'new_c11_tetramer.pdb' == infile:
        '''The C11 nucleosome tetramer with all the protein tails'''
        dna_segnames = ['DNA1', 'DNA2']
        dna_resids.append([1, 693])
        dna_resids.append([693, 1])
        # continuous flexible residues on the first DNA strand
        # recall that range(a, b) excludes upper lim: [a, b)
        flex_resids = [range(1, 31), range(167, 198), range(334, 365),
             range(501, 532), range(667, 694)]
        rigid_groups.append(['A0', 'B0', 'C0', 'D0',
                             'E0', 'F0', 'G0', 'H0'])
        rigid_groups.append(['A1', 'B1', 'C1', 'D1',
                             'E1', 'F1', 'G1', 'H1'])
        rigid_groups.append(['M1', 'N1', 'O1', 'P1',
                             'Q1', 'R1', 'S1', 'T1'])
        rigid_groups.append(['M0', 'N0', 'O0', 'P0',
                             'Q0', 'R0', 'S0', 'T0'])
    elif 'new_dsDNA.pdb' == infile:
        # linker dna file
        try:
            dna_segnames = variables['dna_segnames'][0]
            dna_resids   = variables['dna_resids'][0]
            flex_resids  = variables['flex_resids'][0]
        except:
            dna_segnames = ['DNA1', 'DNA2']
            dna_resids.append([1, 30]) # DNA base pairing
            dna_resids.append([30, 1]) # DNA base pairing
            flex_resids = [range(5, 26)]
    else:
        try:
            dna_segnames = variables['dna_segnames'][0]
            dna_resids   = variables['dna_resids'][0]
            flex_resids  = variables['flex_resids'][0]
            rigid_groups = variables['rigid_groups'][0]
        except:  
            message = "\n~~~ ERROR, unknow pdb file input ~~~\n"
            print_failure(message, txtOutput) 
            return

    i_loop = 0
    remaining_trials = trials
    aa_dcdfiles = []
    cg_dna_dcdfiles = []    
    rigid_dcdfiles = []  
    dna_bead_masks = []
    aa_dna_mask = []
    pkl_file = False
    simulation = []
    while remaining_trials > 0:

        tic = time.time()
        (cg_dna, aa_dna, rigid_mol, vecXYZ, trialbeads, beadgroups, 
         rigid_move_masks, all_beads, dna_bead_masks, rigid_group_masks, 
         rigid_group_mols, aa_all, rigid_mask, aa_dna_mask, bp_per_bead, 
         pkl_file) = get_cg_parameters(flex_resids, dna_resids, dna_segnames, 
                                       infile, refpdb, bp_per_bead, txtOutput, 
                                       rigid_groups, path, rm_pkl, debug, 
                                       i_loop, pkl_file)
        toc = time.time() - tic 
        print 'Total coarse-grain time = %0.3f seconds' % toc        
            
        loop_trials = 100
        if remaining_trials < loop_trials:
            loop_trials = remaining_trials
        remaining_trials -= loop_trials  # increment for the number of trials
        
        if debug:
            print 'loop_trials =', loop_trials
        
        tic = time.time()     
        aa_ofile, cg_dna_ofile, rigid_ofile = dna_mc(loop_trials, i_loop,
            theta_max, theta_z_max, debug, goback, n_dcd_write, keep_unique, 
            keep_cg_files, softrotation, write_flex, runname, ofile, cg_dna, 
            aa_dna, rigid_mol, vecXYZ, lp, trialbeads, beadgroups, 
            rigid_move_masks, all_beads, dna_bead_masks, rigid_group_masks,
            rigid_group_mols, aa_all, rigid_mask, aa_dna_mask, seed)
        aa_dcdfiles.append(aa_ofile)
        cg_dna_dcdfiles.append(cg_dna_ofile)
        rigid_dcdfiles.append(rigid_ofile)
        toc = time.time() - tic
        print 'loop time = %0.3f seconds' % toc
        
        if remaining_trials > 0:
            tic = time.time()
            infile, simulation = minimize(aa_ofile, refpdb, path, psffile,
                                          runname, temperature, openmm_min, 
                                          debug, i_loop, simulation)
            toc = time.time() - tic
            print 'minimization time = %0.3f secords' % toc
            i_loop += 1
            
    # combine the output dcd's
    combine_output(runname, refpdb, txtOutput, ofile, aa_dcdfiles, debug, 
                   keep_cg_files, cg_dna_dcdfiles, rigid_dcdfiles, 
                   rigid_move_masks)
        
    print '\nFinished %d successful DNA MC moves! \n\m/ >.< \m/' % trials

def combine_output(runname, refpdb, txtOutput, ofile, aa_dcdfiles, debug, 
                   keep_cg_files, cg_dna_dcdfiles, rigid_dcdfiles, 
                   rigid_move_masks):
    '''
    this method serves to combine the un-minimized output structures into one 
    complete dcd file
    '''
    output_log_file = open(runname+'/dna_mc/combine_output.log', 'w')
    output_path = runname + '/'
    merge_dcd_files(ofile[:-4], refpdb, 
                    aa_dcdfiles, output_path, output_log_file, txtOutput)
    if keep_cg_files:
        merge_dcd_files('cg_dna', cg_dna_dcdfiles[0][:-4] + '.pdb', 
                    cg_dna_dcdfiles, output_path, output_log_file, txtOutput)
        if len(rigid_move_masks) > 0:       #check that there are actually proteins
            merge_dcd_files('cg_pro', rigid_dcdfiles[0][:-4] + '.pdb', 
                    rigid_dcdfiles, output_path, output_log_file, txtOutput)
    if not debug:
        for aa_ofile in aa_dcdfiles:
            try:
                os.system('rm ' + aa_ofile)
            except:
                message = '\nfailed to remove %s\n' % aa_ofile
                print message
                txtOutput.put(message)
        if keep_cg_files:
            cg_dcdfiles = cg_dna_dcdfiles + rigid_dcdfiles
            for cg_ofile in cg_dcdfiles:
                try:
                    os.system('rm ' + cg_ofile)
                except:
                    message = '\nfailed to remove %s\n' % cg_ofile
                    print message
                    txtOutput.put(message)
        
    outdir = runname + '/dna_mc/'
    try:
        os.system('mv ' + output_path + 'generate/* ' + outdir)
        output_log_file.write('moved all output files from %s to %s\n' % (
            output_path+'generate/', outdir))
        os.system('rmdir ' + output_path + 'generate/')
    except:
        message = ('cannot move combined output dcd files into the project '
                   'directory: %s' % outdir)
        print_failure(message, txtOutput)

    output_log_file.close()

class my_variables(object):
    def __init__(self, pdbfile=None, psffile=None, topfile=None, parmfile=None, 
                 integrator=None):
        pass
        

def minimize(aa_dcd, refpdb, path, psffile, runname, temperature, openmm_min,
             debug, i_loop, simulation):
    '''
    This method is passes the coordinates for the last structure to the open-mm
    minimizer then return the minimized structure
    '''

    max_steps = 5000
    number_of_equilibration_steps = 1000
    number_of_nvt_steps = 1000
    
    toppath = sasconfig._bin_path+'/toppar/'
    testpath = '/Users/curtisj/Desktop/august_sassie_development/svn_utk/sassie_1.0/trunk/testing/'
    testpath = '/home/curtisj/svn_utk/svn/sassie_1.0/trunk/testing/'
    
    variables = my_variables()
    variables.pdbfile = path + refpdb

    # load the pdb then find and load the last structure from the dcd
    raw_structure = sasmol.SasMol(0)
    raw_structure.read_pdb(variables.pdbfile)
    tmp_mol = sasmol.SasMol(0)
    dcdfile = tmp_mol.open_dcd_read(aa_dcd)
    last_frame = dcdfile[2]
    raw_structure.read_single_dcd_step(aa_dcd, last_frame)

    if openmm_min:
        # energy_convergence = None# energy_convergence = None
        energy_convergence = 1.0*unit.kilojoule/unit.mole  # the default value
        if i_loop == 0:
            variables.psffile = psffile
            variables.topfile = toppath+'top_all27_prot_na.inp'
            variables.parmfile = toppath+ 'par_all27_prot_na.inp'            
            step_size = 0.002*unit.picoseconds  # 1000000*0.002 = 2 ns total simulation time
            variables.integrator = LangevinIntegrator(temperature*unit.kelvin, 1/unit.picosecond, step_size)
            simulation = omm.initialize_system(variables)   
        omm.energy_minimization(simulation,raw_structure,energy_convergence,max_steps)
    else:
        print 'no minimization selected, returing raw file'
        return aa_dcd, simulation
    
    min_file = '%s_min.dcd' % aa_dcd[:-4]
    raw_structure.write_dcd(min_file)
    print 'minimized structure: %s' % min_file


    return min_file, simulation

def get_cg_parameters(flex_resids, dna_resids, dna_segnames, infile, refpdb, 
                      bp_per_bead, txtOutput, rigid_groups=[], path='./', rm_pkl=False, 
                      debug=False, i_loop=0, pkl_file=False):
    '''
    This method is designed to generate the coarse-grained run parameters then
    save them to a pickle file for future use.  If the pickle file already 
    exists, it will load that file then compare the new input to see if all 
    parameters meet the current input.  If there are critical differences, it 
    will re-generate the parameters and update the pickle file.  This process 
    drastically improves the load time when running the same setup many times.
    To force the program to regenerate the parameters, set 
    variables['rm_pkl'][0] = True
    to force the removal of the pickle file.  This should occur when 
    re-coarse-graining after a minimization.  
    '''
    try:
        bp_per_bead_val0 = bp_per_bead[0]
        bp_per_bead_array = bp_per_bead
    except:
        bp_per_bead_val0 = bp_per_bead
        
    if not pkl_file:
        pkl_file = path + infile[:-3] + 'pkl'

    if os.path.isfile(pkl_file) and rm_pkl and i_loop == 0:
        print '>>> removing cg paramaters for %s: %s' % (infile, pkl_file)
        os.remove(pkl_file)

    do_rigid = do_cg_dna = do_dna_flex = load_infile = False

    # if pickle_file exists:
    if os.path.isfile(pkl_file):
        print 'loading cg parameters for %s from: %s' % (infile, pkl_file)
        # load pickle_file
        pkl_in = open(pkl_file, 'rb')

        # coordinate dependent
        cg_dna            = pickle.load(pkl_in)
        aa_dna            = pickle.load(pkl_in)
        rigid_mol         = pickle.load(pkl_in)
        aa_all            = pickle.load(pkl_in)
        vecXYZ            = pickle.load(pkl_in)

        # independent of coordinates
        trialbeads        = pickle.load(pkl_in)
        beadgroups        = pickle.load(pkl_in)

        rigid_move_masks  = pickle.load(pkl_in)
        all_beads         = pickle.load(pkl_in)
        dna_bead_masks    = pickle.load(pkl_in)
        rigid_group_masks = pickle.load(pkl_in)
        rigid_group_mols  = pickle.load(pkl_in)
        rigid_mask        = pickle.load(pkl_in)
        aa_dna_mask       = pickle.load(pkl_in)
        
        # input parametrs used to generate these cg-parameters
        infile_old        = pickle.load(pkl_in)
        refpdb_old        = pickle.load(pkl_in)

        dna_resids        = pickle.load(pkl_in)
        dna_segnames      = pickle.load(pkl_in)
        flex_resids_old   = pickle.load(pkl_in)
        rigid_groups_old  = pickle.load(pkl_in)
        bp_per_bead_old   = pickle.load(pkl_in)
        pkl_in.close()
        
        # check if input parameters have changes since last using this pdb
        if rigid_groups != rigid_groups_old or i_loop > 0:
            do_rigid = True
            print '>>>Re-coarse-graining the Protein'               

        if bp_per_bead_val0 != bp_per_bead_old[0] or i_loop > 0:
            do_cg_dna = True
            print '>>>Re-coarse-graining the DNA'       
        else:
            bp_per_bead_array = bp_per_bead_old

        if flex_resids != flex_resids_old:
            do_dna_flex = True
            print '>>>Re-identifying flexible residues'       
            
    else:
        do_cg_dna = do_rigid = load_infile = True
        dna_bead_masks = []
        aa_dna_mask = []
        
    if i_loop > 0 or load_infile:
        # load in the all atom pdb
        aa_all = sasmol.SasMol(0)
        if infile[-3:] == 'dcd':
            aa_all.read_pdb(path + refpdb)  
            tmp_mol = sasmol.SasMol(0)
            dcdfile = tmp_mol.open_dcd_read(infile)
            last_frame = dcdfile[2]
                
            if last_frame <= 0:
                message = 'input dcd file contains no frames or is corrupt\n'
                message += 'using pdb file'
                txtOutput.put(message)
            else:
                aa_all.read_single_dcd_step(infile, last_frame)
            
        elif infile[-3:] == 'pdb':
            aa_all.read_pdb(path + infile)  
            
        else:
            'unknown file type for input structure'
            ruturn

    if do_cg_dna:
        #~~~ DNA ONLY SECTION ~~~#
        (aa_dna, aa_dna_mask, cg_dna, all_beads, dna_bead_masks, 
         vecXYZ, bp_per_bead_array) = make_cg_dna(dna_segnames, dna_resids, 
                                            bp_per_bead_val0, aa_all, txtOutput,
                                            dna_bead_masks, aa_dna_mask)
        do_dna_flex = True
        
    if do_dna_flex:
        #~~~ Determine flexible DNA beads from user input ~~~#
        beadgroups, trialbeads = is_bead_flexible(flex_resids, cg_dna.natoms(), 
                                    dna_resids[0], bp_per_bead_array, debug)

    if do_rigid:
        #~~~ Non-DNA section ~~~#
        (rigid_mol, rigid_mask, rigid_group_masks, rigid_group_mols, rigid_move_masks
         ) = make_rigid_groups(aa_all, rigid_groups)
        
    if do_cg_dna or do_rigid or do_dna_flex:
        print 'cg %s using updated parameters. Result saved to: %s' % (infile,
                                                                       pkl_file)

        # create the pickle_file for future use
        pkl_out = open(pkl_file, 'wb')

        # coordinate dependent
        pickle.dump(cg_dna, pkl_out, -1)
        pickle.dump(aa_dna, pkl_out, -1)
        pickle.dump(rigid_mol, pkl_out, -1)
        pickle.dump(aa_all, pkl_out, -1)
        pickle.dump(vecXYZ, pkl_out, -1)

        # independent of coordinates
        pickle.dump(trialbeads, pkl_out, -1)
        pickle.dump(beadgroups, pkl_out, -1) 

        pickle.dump(rigid_move_masks, pkl_out, -1) 
        pickle.dump(all_beads, pkl_out, -1) 
        pickle.dump(dna_bead_masks, pkl_out, -1)
        pickle.dump(rigid_group_masks, pkl_out, -1)  
        pickle.dump(rigid_group_mols, pkl_out, -1)  
        pickle.dump(rigid_mask, pkl_out, -1)  
        pickle.dump(aa_dna_mask, pkl_out, -1)    

        # input parameters used to get these cg-parameters
        pickle.dump(infile, pkl_out, -1)        
        pickle.dump(refpdb, pkl_out, -1)        
        
        pickle.dump(dna_resids, pkl_out, -1)
        pickle.dump(dna_segnames, pkl_out, -1)        
        pickle.dump(flex_resids, pkl_out, -1)
        pickle.dump(rigid_groups, pkl_out, -1)        
        pickle.dump(bp_per_bead_array, pkl_out, -1)
        pkl_out.close()
    
    # this is important for re-aligning the proteins after moving them
    # cg_pro_orig = sasmol.SasMol(0)
    # error = cg_pro.copy_molecule_using_mask(cg_pro_orig,
                                        # numpy.ones(len(cg_pro.coor()[0])), 0)
    
    return (cg_dna, aa_dna, rigid_mol, vecXYZ, trialbeads, beadgroups, 
            rigid_move_masks, all_beads, dna_bead_masks, rigid_group_masks,
            rigid_group_mols, aa_all, rigid_mask, aa_dna_mask, 
            bp_per_bead_array, pkl_file)

def prepare_dna_mc_input(variables):
    # ~~~~Generate 'flex_resid' list~~~~
    theta_max             = variables['theta_max'][0]
    theta_z_max           = variables['theta_z_max'][0]
    n_flex_regions        = variables['n_flex_regions'][0]
    first_res_per_region  = variables['first_res_per_region'][0]
    n_cont_res_per_region = variables['n_cont_res_per_region'][0]

    assert len(theta_max) == n_flex_regions, 'theta_max should have %d values' % n_flex_regions
    assert len(theta_z_max) == n_flex_regions, 'theta_z_max should have %d values' % n_flex_regions
    assert len(first_res_per_region) == n_flex_regions, 'first_res_per_region should have %d values' % n_flex_regions
    assert len(n_cont_res_per_region) == n_flex_regions, 'n_cont_res_per_region should have %d values' % n_flex_regions
    
    flex_resids = []
    for i in xrange(n_flex_regions):
        first = first_res_per_region[i]
        last  = first_res_per_region[i] + n_cont_res_per_region[i]
        flex_resids.append(range(first, last))

    variables['flex_resids'] = (flex_resids, 'list_of_lists')

    # ~~~~Generate 'pro_groups' list~~~~
    n_rigid_groups = variables['n_rigid_groups'][0]
    assert n_rigid_groups <= n_flex_regions, '# of protein groups must be <= to # of flexible DNA regions'
    rigid_groups = []
    for i in xrange(n_rigid_groups):
        rigid_groups.append(variables['rigid_group'+str(i+1)][0])
    variables['rigid_groups'] = (rigid_groups, 'list_of_lists')
        
    # ~~~~Generate 'dna_resids' list~~~~
    dna_resids = [variables['dna1_resids'][0], variables['dna2_resids'][0]]
    variables['dna_resids'] = (dna_resids, 'list_of_lists')
    
    return variables


if __name__ == "__main__":

    ''' wanted to implement this for error handling but never did'''

    s_theta_max = '20'
    if 'theta_z_max' in locals():
        s_theta_z_max = ''
        for (i, theta) in enumerate(theta_z_max):
            if i > 0:
                s_theta_z_max += ', '
            s_theta_z_max += str(theta)    
    else:
        s_theta_z_max = s_theta_max
        
    svariables = {}
    
    # User Input
    svariables['runname'] = ('testing', 'string')
    svariables['path']    = ('./', 'string')
    svariables['infile']  = ('new_dsDNA60.pdb', 'string')
    svariables['refpdb']  = ('new_dsDNA60.pdb', 'string')
    svariables['psffile'] = ('new_dsDNA60.psf', 'string')
    svariables['ofile']   = ('new_dsDNA60_mc.dcd', 'string')
    svariables['trials']  = ('25', 'int')
    svariables['goback']  = ('50', 'int')
    
    # Molecule Specific Input
    svariables['n_flex_regions']        = ('1', 'int')
    svariables['theta_max']             = ('15', 'float_array')
    svariables['theta_z_max']           = ('5', 'float_array') # provide 's_theta_z_max' to use a different theta max for twisting vs bending
    svariables['first_res_per_region']  = ('16', 'int_array')
    svariables['n_cont_res_per_region'] = ('30', 'int_array')
    svariables['align_low_res']         = ('1', 'int')  # not yet implemented
    svariables['align_high_res']        = ('15', 'int')  # not yet implemented
    svariables['dna_segnames']          = ('DNA1, DNA2', 'string_array')  
    svariables['dna1_resids']           = ('1, 60', 'int_array')
    svariables['dna2_resids']           = ('120, 61', 'int_array')
    svariables['n_pro_groups']          = ('0', 'int')
    svariables['pro_group1']            = ('A0, B0, C0, D0, E0, F0, G0, H0', 'string_array')
    svariables['pro_group2']            = ('', 'string_array') # place holder
    
    # Specialized/Advanced Inputs
    svariables['bp_per_bead']   = ('1', 'int') # set to N > 1 to have N base-pairs coarse-grained into 1 bead
    svariables['softrotation']  = ('1', 'int') # set to N > 1 to apply rotations averaged over N coars-grained beads
    svariables['n_dcd_write']   = ('1', 'int') # set to N > 1 to only record structures after every N trials
    svariables['seed']          = ('0', 'int') # goback random seed
    svariables['temperature']   = ('300', 'float') # may be incorperated for electrostatics
    svariables['debug']         = ('False', 'bool') # 'True' will display debug output
    svariables['write_flex']    = ('False', 'bool') # 'True' will generate a file labeling paird DNA base pairs for all flexible beads (useful for calculating dihedral angles)
    svariables['keep_cg_files'] = ('True', 'bool') # 'True' will keep the coarse-grain DNA and protein pdb and dcd files
    svariables['keep_unique']   = ('True', 'bool') # 'False' will keep duplicate structure when move fails
    svariables['rm_pkl']        = ('False', 'bool') # 'True' will remove the coarse-grain pkl file forcing a re-coarse-graining (can be time consuming often necessary)
    svariables['openmm_min']    = ('False', 'bool') # 'True' will remove the coarse-grain pkl file forcing a re-coarse-graining (can be time consuming often necessary)    

    error, variables = input_filter.type_check_and_convert(svariables)
    variables = prepare_dna_mc_input(variables)

    main(variables)
