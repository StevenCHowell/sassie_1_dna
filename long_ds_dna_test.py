#!/usr/bin/python
# $Id$
# time using FORTRAN double loop, N=1000, iters=1000 (so 1*10^6 steps): 958.887075186 seconds
# time using python double loop, N=1000, iters=1000 (so 1*10^6 steps): 

import sassie.sasmol.sasmol as sasmol
import numpy as np, string, os, locale, sys, random, time

def write_xyz(filename,coor,comment,frame):

    natoms=len(coor)
    if(frame==0):
        outfile=open(filename,'w')
    else:
        outfile=open(filename,'a')

    outfile.write("%i\n" % natoms)
    outfile.write('%s\t%s\n' % ('B',str(comment)))
    for i in range(natoms):
        outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
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

    # print 'natoms = ',natoms

    cg_dna = sasmol.SasMol(0)
    basis_filter = "(resname[i] == 'ADE' or resname[i] == 'GUA') and (name[i] == 'N1' and (resid[i] < 21) and segname[i]=='DNA1')"
    error,mask = aa_dna.get_subset_mask(basis_filter) 

    frame = 0
    error = aa_dna.copy_molecule_using_mask(cg_dna,mask,frame) 
    error,coor=aa_dna.get_coor_using_mask(frame,mask)

    cg_dna.setCoor(coor)

    diameter = 6.0*3.4

    cg_natoms = cg_dna.natoms()
#        print 'cg_natoms = ',cg_natoms
    new_coor = np.zeros((1,cg_natoms,3),np.float)
#        print new_coor

    for i in xrange(cg_natoms):
#                print "i = ",i
        new_coor[0][i][2] = i*diameter

    cg_dna.setCoor(new_coor)
#        print new_coor
#        print cg_dna.coor()

#        charge = cg_dna.beta()
#        cg_dna.setCharge(charge)

#        print 'loc = ',cg_dna.loc()
#        print 'rescode = ',cg_dna.rescode()
#        print 'charge = ',cg_dna.charge()

#        print 'type(index) = ',type(cg_dna.index())
    index = cg_dna.index()

#sh        coor = cg_dna.coor()[0]
#sh        comment = 'test'
#sh        write_xyz('cg_dna.xyz',coor,comment,frame)

    infile=open('dum.pdb','w')
    for i in xrange(cg_dna.natoms()):
        this_index = index[i]

        sx = cg_dna._coor[frame,i,0]#[:8]
        sy = cg_dna._coor[frame,i,1]#[:8]
        sz = cg_dna._coor[frame,i,2]#[:8]

        infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s\n" % (cg_dna.atom()[i],this_index,cg_dna.name()[i],cg_dna.loc()[i],cg_dna.resname()[i],cg_dna.chain()[i],cg_dna.resid()[i],cg_dna.rescode()[i],sx,sy,sz,cg_dna.occupancy()[i],cg_dna.beta()[i],cg_dna.segname()[i],cg_dna.element()[i],cg_dna.charge()[i]))

    infile.write('END\n')
    infile.close()
    os.remove('dum.pdb') # remove the temporary pdb file

    cg_dna.write_pdb("cg_test.pdb",frame,'w')

#sh        cg_dna.write_dcd("cg_test.dcd")

    return cg_dna

def make_cgDNA_model(all_atom_pdb, chain1, chain2, resid1, resid2, flexResid, bp_perBead):
    '''
    this function creates a cg bead model from DNA chains supplied as input
    make_cgDNA_model takes specific DNA to coarse grain as input
    '''
    aa_dna = sasmol.SasMol(0)
    aa_dna.read_pdb(all_atom_pdb)
    natoms = aa_dna.natoms() #; print 'natoms = ',natoms

    cg_dna = sasmol.SasMol(0)  # need to give this propertios

    # check the input
    #        assert resid1.size == resid2.size == 2*chain1.size == 2*chain2.size == 2*bp_perBead.size, "missing an input value for coarse-graining the DNA"
    assert np.abs(resid1[1]-resid1[0]) == np.abs(resid2[1]-resid2[0]), "number of paired bases in chains are not the same"

    frame = 0
    bps = resid1[1]-resid1[0]+1
    nbeads = int(np.round(bps/float(bp_perBead))) #; print 'nbeads = ',nbeads

    # populate the cg_dna sasmol object with atom properties
    s = 0
    i = 0
    while s < nbeads:
        basis_filter = "((name[i] == 'P') and (resid[i] < "+str(nbeads+i)+"))" #; print 'bf = ',basis_filter
        error, mask = aa_dna.get_subset_mask(basis_filter)   #; print mask
        s = np.sum(mask)
        i += 1

    error = aa_dna.copy_molecule_using_mask(cg_dna,mask,frame) 

    cg_coor = np.zeros((1,nbeads,3)) # initialize array to store the coordinates of the beads

    # select the cg beads coordinates from the PDB
    #s if resid1[1] < resid1[0]:
    #s         tmp = resid1[i,1]
    #s         resid1[1] = resid1[0]
    #s         resid1[0] = tmp
    #s if resid2[1] < resid2[0]:
    #s         tmp = resid2[1]
    #s         resid1[1] = resid2[0]
    #s         resid2[0] = tmp

    npairs = resid1[1]-resid1[0]+1
    r1a = resid1[0]
    r2a = resid2[0]
    allBeads = []    # list to contain all the all atom sasmol object of each bead 
    flex = np.zeros(nbeads)

    for j in xrange(nbeads):
        if j+1 == nbeads:
            # remain = bps%bp_perBead # not sure if this is the best way to handle the remainder
            bp_perBead = npairs - bp_perBead*j  # to account for the remainder

        bead = sasmol.SasMol(0)
        if resid1[0] < resid1[1]:
            r1b = r1a + bp_perBead
        else:
            r1b = r1a - bp_perBead
        if resid2[0] < resid2[1]:
            r2b = r2a + bp_perBead
        else:
            r2b = r2a - bp_perBead

        basis_filter = "((resid[i] > "+str(r1a-1)+" and resid[i] < "+str(r1b)+") and (chain[i]=='"+chain1+"')) or ((resid[i] > "+str(r2a-1)+" and resid[i] < "+str(r2b)+") and (chain[i]=='"+chain2+"'))"

        error, mask = aa_dna.get_subset_mask(basis_filter)   # create a mask to select the atoms for the bead
        error = aa_dna.copy_molecule_using_mask(bead,mask,frame) 
        com = bead.calccom(frame)
        cg_coor[0,j,:] = com         # store the coordinates of the com for the bead coordinate
        # print 'com ', j, '= ', com

        bead.center                  # calculate the coodinates of the atoms in the bead with reference to the bead com
        allBeads.append(bead)

        # check if this is a flexible bead
        if [i for i in [r1a, r1b-1, r2a, r2b-1] if i in flexResid]:
            flex[j] = 1
            # print 'bead ', j, 'is flexible'
            # print r1a, r1b-1, r2a, r2b-1
            # print flexResid

        # setup for next iteration
        r1a = r1b
        r2a = r2b

    # set the bead coordinates to the calculated com coordinates
    #print 'cg_coor = \n', cg_coor
    flexible = np.zeros(np.sum(flex))
    j = 0
    for i in xrange(flex.size):
        if flex[i]:
            flexible[j] = i
            j+=1

    print 'flex = \n', flex
    print 'flexible = ',flexible



    cg_dna.setCoor(cg_coor)

    cg_dna.write_pdb("cgDNA.pdb",frame,'w')

    vecXYZ = np.zeros((nbeads*3,3))
    vecXYZ[0:nbeads] = [1,0,0]
    vecXYZ[nbeads:2*nbeads] = [0,1,0]
    vecXYZ[2*nbeads:3*nbeads] = [0,0,1]
    return (cg_dna, vecXYZ, allBeads, flexible)

def recover_aaDNA_model(cg_dna, vecXYZ_f, allBeads):
    cg_natoms = cg_dna.natoms()
    coor = cg_dna.coor()[0] #; print 'coor = ', coor

    av = checkMag(vecXYZ) #; print 'avMag = ',av

    # split vecXYZ into three matrices
    vecX = vecXYZ_f[0:cg_natoms]
    vecY = vecXYZ_f[cg_natoms:2*cg_natoms]
    vecZ = vecXYZ_f[2*cg_natoms:3*cg_natoms]

    for i in xrange(cg_natoms):
        # R is the matrix that would align the rotated coordinates back to the original
        R = align2xyz(vecX[i,:], vecY[i,:], vecZ[i,:]) #; print 'R = ', R

        # M will rotate cooridates from the original reference to the rotated coordinates of the bead then translate the com to the beads com
        M = R.transpose() 
        M[3,:3] = coor[i,:] #; print 'M = ', M
        (r,c) = allBeads[i].coor()[0].shape
        beadCoor = np.ones((r,4)) 
        beadCoor[:,:3] = allBeads[i].coor()[0] #; print 'beadCoor = ',beadCoor

        allBeads[i].setCoor(np.dot(beadCoor,M))

    return allBeads

def align2xyz(vecX, vecY, vecZ):

    tmp_coor = np.zeros((2,4))
    tmp_coor[:,3] = 1
    tmp_coor[1,0:3] = vecZ
    A1 = align2z(tmp_coor)

    newX = np.dot(vecX,A1[0:3,0:3]) ;# print 'newX = ', newX
    #newY = np.dot(vecY,A1[0:3,0:3]) ;# print 'newY = ', newY
    assert newX[2] < 10e-5, "ERROR!!! z-component of newX is not zero and it should be"

    thetaZ_x = -np.arctan2(newX[1],newX[0])  ;#  print 'thetaZ_x = ', thetaZ_x
    #thetaZ_y = -np.arctan2(newY[1],-newY[0])  ;#  print 'thetaZ_y = ', thetaZ_y

    A2 = rotate4x4('z', thetaZ_x)  #; print 'A2 = ', A2

    A = np.dot(A1,A2)

    newY = np.dot(vecY,A[0:3,0:3]) ;# print 'finalY = ', newY
    assert newY[0]+newY[2] < 1+10e-5, "ERROR!!! newY is not aligned to the y-axis and it should be"

    return A

def rotate4x4(axis, theta):
    R = np.eye(4)
    ct = np.cos(theta)
    st = np.sin(theta)
    if axis.lower()=='x':
        (R[1,1], R[1,2]) = (ct,  st)
        (R[2,1], R[2,2]) = (-st, ct)
    elif axis.lower()=='y':
        (R[0,0], R[0,2]) = (ct, -st)
        (R[2,0], R[2,2]) = (st,  ct)
    elif axis.lower()=='z':
        (R[0,0], R[0,1]) = ( ct, st)
        (R[1,0], R[1,1]) = (-st, ct)
    else:
        assert True, "ERROR!!! did not recognize rotation axis"

    # print 'R = ', R

    return R


def move2origin(coor4):
    T    = np.eye(4,dtype=np.float)
    Ti   = np.copy(T)
    #s print 'coor passed to T, Ti:\n',coor4
    #s print 'size = ', coor4.shape
    if not coor4.shape < (1,4):
        T[3,0:3] = -coor4[0,0:3]
        Ti[3,0:3] = coor4[0,0:3]

    #s print 'T:\n',T
    #s print 'test:\n',np.dot(coor4,T)
    #s print 'Ti:\n',Ti
    #s print 'testI:\n',np.dot(np.dot(coor4,T),Ti)
    return (T, Ti)

def align2z(coor4):
    '''
    function designed to align the axis connecting the first 2 coordinates of an array 
    of (1,4) coodinate vectors to the z-axis
    '''
    A = np.eye(4,dtype=np.float)
    #s Axz  = np.eye(4,dtype=np.float)
    #s Az   = np.eye(4,dtype=np.float)
    assert all(coor4[0] == [0., 0., 0., 1.,]), "coordinates passed to align2z were not translated to the origin"

    if coor4.shape > (1, 4):
        (u, v, w) = coor4[1,0:3]
        small = 1E-14 # small limit to make sure not dividing by zero

        d1 = np.sqrt(u**2+v**2)
        d2 = np.sqrt(u**2+v**2+w**2)
        if d1 > small:
            (A[0][0], A[0][1], A[0][2]) = (u/d1*w/d2, -v/d1, u/d2)
            (A[1][0], A[1][1], A[1][2]) = (v/d1*w/d2, u/d1, v/d2)
            (A[2][0], A[2][1], A[2][2]) = (-d1/d2, 0, w/d2)   #; print 'A1:\n',A

        #s # align to the x-z plane
        #s if v > small:           # check if already aligned to xz-plane
        #s         (Axz[0][0], Axz[0][1]) = (u/d1, -v/d1)
        #s         (Axz[1][0], Axz[1][1]) = (v/d1,  u/d1)  #; print 'Axz= \n', Axz
        #s 
        #s # align to the z-axis
        #s if d1 > small:
        #s         (Az[0][0], Az[0][2]) = (w/d2,  d1/d2)
        #s         (Az[2][0], Az[2][2]) = (-d1/d2, w/d2)  #; print 'Az= \n', Az
    else:
        print 'no point to align'
    #s A = np.dot(Axz,Az) #;   print 'A3=\n',A
    return A


def beadRotate(coor3,vecX,vecY,vecZ,thetas,nSoft):
    '''
    this function is designed to generate a modified version of the input coordinates (coor3)
    it moves the all coordinates so the first is at the origin
    aligns the first two coordinates to be along the z-axis
    then rotates all coordinates successivly by
       1- thetas[2] about the z-axis
       2- thetas[0] about the x-axis
       3- thetas[1] about the y-axis
    it then undoes the alignment to the z-axis and translates all coordinates so the first coordinate
    is where it started
    '''
    # add a fourth column of ones to each bead origin and orientation vectors (this incorporates translations)
    (natoms,col) = coor3.shape

    # make sure the number of cg beads and orientation vectors match (this will not include any rigid components)
    assert (natoms,col) == vecX.shape == vecY.shape == vecZ.shape, "different number of bead origins and orientations"

    # initialize vector arrays for coordinates and orientation vectors
    coor4 = np.ones((natoms,4),np.float)
    X = np.copy(coor4)
    Y = np.copy(coor4)
    Z = np.copy(coor4)
    coor4[:,0:3] = coor3     #; print 'coor4 = ',coor4
    X[:,0:3] = vecX
    Y[:,0:3] = vecY
    Z[:,0:3] = vecZ

    # create the translation-rotation matrix
    # This is intended to be multiplied from the right (unlike standard matrix multiplication) so as not to require transposing the coordinate vectors.
    # print '(tx, ty, tz) in degrees:', thetas
    tx = thetas[0]*(np.pi/180.0)
    ty = thetas[1]*(np.pi/180.0)
    tz = thetas[2]*(np.pi/180.0)
    # print '(tx, ty, tz) in radians:',(tx,ty,tz)
    cx = np.cos(tx)
    sx = np.sin(tx)
    cy = np.cos(ty)
    sy = np.sin(ty)
    cz = np.cos(tz)
    sz = np.sin(tz)

    # initialize the transformation pieces
    # consolidated method of defining the rotation matrices
    R = np.eye(4,dtype=np.float)
    (R[0][0], R[0][1], R[0][2]) = ( cy*cz, cy*sz, -sy )
    (R[1][0], R[1][1], R[1][2]) = ( sx*sy*cz-cx*sz, sx*sy*sz+cx*cz, sx*cy )
    (R[2][0], R[2][1], R[2][2]) = ( sx*sz+cx*sy*cz, cx*sy*sz-sx*cz, cx*cy )  #; print 'R:\n', R

    # verbose method of defining the rotation matrices (slower)
    #s Rx = np.eye(4,dtype=np.float)
    #s Ry = np.eye(4,dtype=np.float)
    #s Rz = np.eye(4,dtype=np.float)
    #s 
    #s (Rx[1][1], Rx[1][2]) = ( cx, sx)
    #s (Rx[2][1], Rx[2][2]) = (-sx, cx)  #; print 'Rx:\n', Rx
    #s 
    #s (Ry[0][0], Ry[0][2]) = (cy, -sy)
    #s (Ry[2][0], Ry[2][2]) = (sy,  cy)  #; print 'Ry:\n', Ry
    #s 
    #s (Rz[0][0], Rz[0][1]) = ( cz, sz)
    #s (Rz[1][0], Rz[1][1]) = (-sz, cz)  #; print 'Rz:\n', Rz
    #s 
    #s R = np.dot(np.dot(Rx, Ry), Rz) #; print "Rxyz = \n", Rxyz

    # print 'original coor:\n', coor4
    (T0, Ti0) = move2origin(coor4)
    coor4 = np.dot(coor4,T0)  #; print 'moved to origin coor:\n', coor4

    A = align2z(coor4)

    coor4 = np.dot(coor4,A)         #; print 'aligned coor to z-axis:\n', coor4
    coor4 = np.dot(coor4,R) #; print 'step 0 rotated about first angle coor:\n', coor4

    # repeat rotation for softening the bend
    for i in xrange(1,nSoft):
        (T, Ti) = move2origin(coor4[i:])
        coor4[i:] = np.dot(coor4[i:],T)  # move to origin
        # print 'moved2origin coor:\n',coor4

        coor4[i:] = np.dot(coor4[i:],R)
        # print 'step %d' %i,'rotated coor:\n',coor4

        coor4[i:] = np.dot(coor4[i:],Ti) # return to original position
        # print 'returned from origin coor:\n',coor4
        # the coarse grained beads local coordinates should not be translated, only rotated
        X[i:] = np.dot(X[i:],R)
        Y[i:] = np.dot(Y[i:],R)
        Z[i:] = np.dot(Z[i:],R)

    coor4 = np.dot(coor4,A.transpose())
    # print 'un-aligned:\n',coor4
    coor4 = np.dot(coor4,Ti0)
    # print 'returned from origin coor:\n',coor4
    # this returns the modified positions and orientations for all but the first (reference) bead
    return (coor4[1:,0:3], X[1:,0:3], Y[1:,0:3], Z[1:,0:3])

def checkU(coor):
    '''
    module to generate and check the u-vectors pointing from bead to bead
    makes sure they are all the same length
    returns the u-vectors and their magnitude (l)
    '''
    erLim = 1e-3
    (r,c) = coor.shape
    u = coor[1:,:]-coor[:r-1,:]                 # u-vectors pointing from bead to bead
    lu = np.sqrt(u[:,0]**2+u[:,1]**2+u[:,2]**2) # magnitues of u-vectors -> distance btwn beads
    l = np.mean(lu)                          # average distance btwn bead 

    #print '\n\ncoor:\n', coor
    #print 'u-vectors:\n', u
    #print 'magnitude of u-vectors:\n', lu
    #print 'average magnitude of u =', l
    #print '\n\n'

    test = np.abs(lu-l) > erLim  # check if any of the l-lengths are different
    #s if test.any():
    #s         print 'ERROR: the beads are not uniformly spaced'
    #s         print 'u = \n',u
    #s         print test
    #s         # print 'difference = ', u[:,2] - l
    #s         print 'distance = ', lu

    return (u, l)

def checkMag(vec):
    '''
    module to check that the magnitude of all vectors is the same
    '''
    erLim = 1e-3
    (r,c) = vec.shape
    m = np.zeros(r)
    for i in xrange(r):
        m[i] = mag(vec[i])

    avMag = np.mean(m)
    test = np.abs(m-avMag) > erLim  # check if any of the l-lengths are different
    if test.any():
        print 'ERROR: the vectors do not have uniformly magnitude'
        print 'm = \n',m
        print test
        # print 'difference = ', u[:,2] - l
        print 'vec = ', vec

    return (avMag)

def energyBend(lpl,u,l):
    ''' 
    this function finds the E/kT for N beads connected by N-1 rods
    lpl is the persistence length divided by the rod length: lp/l
    u contains the N-1 unit vectors representing each rod
    '''

    (nu,col) = u.shape          # nu = nbeads - 1
    uk0 = u[:nu-1,:]/l          #; print 'u_k= \n',uk0
    uk1 = u[1:,:]/l             #; print 'u_k= \n',uk1
    #s print 'uk0*uk1\n:',uk0*uk1

    # checkMag(u)       # for debugging

    res = (nu-1) - np.sum(uk0*uk1)  # calculate the sum of their products

    return res*lpl

def energyWCA(w,coor,wca0,trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    '''
    wca1 = np.copy(wca0)
    test = 2.**(1./6.)*w
    (N, col) = coor.shape
    for i in xrange(trial_bead,N):
        ri = coor[i,:]
        for j in xrange(0,i):
            rj = coor[j,:]
            rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)
            #s print '(2^(1/6)*w, rij) = ', (test, rij)
            if rij < test:
                wca1[i,j] = (w/rij)**12.-(w/rij)**6.+0.25
    res = 4*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def FenergyWCA(w,coor,wca0,trial_bead):
    '''
    this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
    beads connected by N-1 rods together representing a worm-like chain.
    w is the effective width of the chain
    coor contains the xyz-coordinates of the beads
    and does all this using FORTRAN
    '''
    import sys ; sys.path.append('./')
    import electrostatics

    wca1 = np.copy(wca0)
    test = 2.**(1./6.)*w
    (N, col) = coor.shape

    wca1 = electrostatics.calcwca(coor,trial_bead+1,w,wca1)  
    #  increment trial bead by one because it will index from 1 instead of 0

    res = 4.*np.sum(wca1)
    #s print 'U_wca =', res*4
    return (res, wca1)

def dna_mc(nsteps,cg_dna,vecXYZ,lp,w,theta_max,flexible,nSoft=3,f=True):
    '''
    this function perform nsteps Monte-Carlo moves on the cg_dna
    '''
    dcdOutFile = cg_dna.open_dcd_write("cg_dna_moves.dcd")        

    nbeads = cg_dna.natoms()  # need to change this so there could be components that are not the beads
    nflex = flexible.size
    coor = np.copy(cg_dna.coor()[0])

    (u, l) = checkU(coor) # get the vectors pointing from bead to bead and make sure they are equidistant

    lpl = lp/l  # setup the presistence length paramater

    # calculate the energy of the starting positions
    wca0 = np.zeros((nbeads,nbeads))
    Ub0 = energyBend(lpl,u,l)
    if f:
        (Uwca0, wca0) = FenergyWCA(w,coor,wca0,0)
    else:
        (Uwca0, wca0) = energyWCA(w,coor,wca0,0)

    U_T0 = Ub0 + Uwca0
    wca1 = np.copy(wca0)

    # split vecXYZ into three matrices
    cg_natoms = cg_dna.natoms()
    vecX = vecXYZ[0:cg_natoms]
    vecY = vecXYZ[cg_natoms:2*cg_natoms]
    vecZ = vecXYZ[2*cg_natoms:3*cg_natoms]

    a = 0 # times configuration was accepted
    r = 0 # times configuration was rejected

    for i in xrange(nsteps):

        trial_bead = flexible(int((nflex-1)*random.random())+1) ## CHECK THIS ##
        # print 'trial_bead =', trial_bead

        thetaz_max = np.float(theta_max)/10. # this should be something small 
        thetaX = theta_max * random.random() - theta_max/2
        thetaY = theta_max * random.random() - theta_max/2
        thetaZ = thetaz_max *random.random() - thetaz_max/2 
        thetaXYZ = [thetaX/nSoft, thetaY/nSoft, thetaZ/nSoft]

        # generate a newly rotated model

        #s for i in xrange(3):
        #s         thetas = np.zeros(3)
        #s         thetas[i] = thetaXYZ[i]
        #s         print 'new thetas:', thetas

        #s print 'coor before:\n',coor
        (coor[trial_bead:],vecX[trial_bead:],vecY[trial_bead:],vecZ[trial_bead:]) = beadRotate(coor[trial_bead-1:],vecX[trial_bead-1:],vecY[trial_bead-1:],vecZ[trial_bead-1:],thetaXYZ,nSoft)
        #s print 'coor after\n',coor

        # calculate the change in energy (dU) and boltzman factor (p) for the new model
        (u, l) = checkU(coor)
        Ub1 = energyBend(lpl,u,l)
        if f:
            (Uwca1,wca1) = FenergyWCA(w,coor,wca0,trial_bead)
        else: 
            (Uwca1,wca1) = energyWCA(w,coor,wca0,trial_bead)

        U_T1 =  Ub1 + Uwca1
        dU = U_T1 - U_T0
        p = np.exp(-dU)
        test = random.random()

        # output the results to a text file
        #s outData.write("%1.3e\t %1.3e\t %1.3e\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %d\t %0.1f\t %0.1f\n" %(dU, Ub1-Ub0, Uwca1-Uwca0, p, test, w, lp, trial_bead, thetaX, thetaY) )
        # print '(dU,dUb,dUwca,p,rand) =', (dU,Ub1-Ub0, Uwca1-Uwca0, p, test)

        # if accepted write new coordinates, else write old again
        #if True:
        if test < p:
            #print 'wrote new dcd frame (end of loop',i,' trial_bead=',trial_bead,')'+' accepted new configuration'
            a += 1
            cg_dna.coor()[0] = np.copy(coor)
            wca0 = np.copy(wca1)
            U_T0 = U_T1
        else :
            #print 'wrote new dcd frame (end of loop',i,' trial_bead=',trial_bead,')'+' rejected new configuration'
            r += 1
            coor = np.copy(cg_dna.coor()[0])   # reset the coordinates

        #s cg_dna.write_dcd_step(dcdOutFile,0,0)

    # close output files
    cg_dna.close_dcd_write(dcdOutFile)
    #s outData.close()

    return (cg_dna, vecXYZ, a, r)

def makeLongDNA(n_lp):
    print 'making DNA that is %d*lp long' %n_lp

    # 15 bp/bead or 51 A/bead (3.4 A/bp)

    lp = 530 # persistence length in A
    l = 2**(1./6.)*46   # separation distance between beads = 51.6A

    longDNA = sasmol.SasMol(0)
    L = n_lp*lp
    N = int(L/l)
    natoms = N+1
    print 'natoms = ',natoms
    longDNA._L = L
    longDNA._natoms = natoms

    longCoor = np.zeros((1,natoms,3),np.float)    # initialize the long DNA coordinates
    longCoor[0][:,2] = range(natoms)      # set the z-values to the index of the array
    longCoor *= l                   # scale the z-values to the right seperation
    # print longCoor[-5:]

    longDNA.setCoor(longCoor)
    longDNA.setElement(['C']*natoms)

    vecXYZ = np.zeros((natoms*3,3))
    vecXYZ[0:natoms] = [1,0,0]
    vecXYZ[natoms:2*natoms] = [0,1,0]
    vecXYZ[2*natoms:3*natoms] = [0,0,1]
    # n = L/l                         # number of times need to repeat the grain
    # print '(l, L, n)', (l, L, n) 

    return (longDNA, vecXYZ)

def mag(vec):
    import numpy as np

    (r,) = vec.shape
    sumSquares = 0.0
    for i in xrange(r):
        sumSquares += vec[i]**2

    return np.sqrt(sumSquares)

def closeAll():
    for x in xrange(100):
        plt.close()

if __name__ == "__main__":

    ### Psuedo code for intended functionality ###
    # Create the bead coordinates
    # Iterate for 10^6 steps (expect to be overkill)
    #   Call ds_dna_method to run for 100,000 trials sample x/y/z separate
    #   Save coordinates to a dcd file
    #   Write the Rg and Re to an output text file
    #   
    #
    #

    # ----- Modify these ---
    iters = 2       # intend to use 1E6
    nsteps = 1000   # intend to use 1E5
    theta_max = np.float(15)
    Llp = 15    # L/lp
    # ----- Modify these ---

    # ----- DO NOT Modify these ---
    lp = 530      # persistence length  (lp = 530A)
    w = 46        # width of chain in A (w = 46A)l
    # ----- DO NOT Modify these ---

    makeLongDNA(Llp)        

    rg0 = cg_dna.calcrg(0)/lp
    re0 = mag(cg_dna.coor()[0,-1]-cg_dna.coor()[0,0])/lp

    if show:
        import pylab
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        # ax = fig.add_subplot(iters, 1, 1, projection='3d')
        ax = Axes3D(fig)
        ax.scatter(cg_dna.coor()[0,:,0],cg_dna.coor()[0,:,1],cg_dna.coor()[0,:,2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fig.suptitle('After %d steps Rg/lp=%f  Re/lp=%f' %(0, rg0, re0))
        plt.show()

    print 'iter:', 0,'of',iters,' (a, r) = (  , )   rg/lp=',rg0, 're/lp=',re0

    timestr = time.strftime("%y%m%d_%H%M%S")
    fName = timestr + '_%dlp_dnaMoves.o' %Llp
    outData = open(fName,'a')
    #iteratons rg/lp a
    outData.write("# L=%d\t iters=%d\t nsteps=%d\t nSoft=%d\t theta_max=%f \n# moves\t rg/lp\t\t re/lp\t\t a\t r\n" %(Llp,iters,nsteps, nSoft, theta_max))
    outData.close()

    tic = time.time()
    for i in xrange(iters):
        (cg_dna,vecXYZ, a, r) = dna_mc(nsteps,cg_dna,vecXYZ,lp,w,theta_max,flexible,nSoft,f)
        rg_lp[i] = cg_dna.calcrg(0)/lp
        re_lp[i] = mag(cg_dna.coor()[0,-1]-cg_dna.coor()[0,0])/lp
        print 'iter:', i+1,'of',iters,' (a, r) = ',(a,r), 'rg/lp=',rg_lp[i], 're/lp=',re_lp[i]

        #output
        outData = open(fName,'a')
        outData.write("%d\t%f\t%f\t%d\t%r\n" %((i+1)*nsteps, rg_lp[i], re_lp[i], a, r) )
        outData.close()
        if show:
            fig = pylab.figure()
            ax = Axes3D(fig)
            ax.scatter(cg_dna.coor()[0,:,0],cg_dna.coor()[0,:,1],cg_dna.coor()[0,:,2])
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            fig.suptitle('After %d steps Rg/lp=%f  Re/lp=%f' %((i+1)*nsteps, rg_lp[i], re_lp[i]))

    toc = time.time() - tic
    print 'run time =',toc,'seconds'

    #recover an all atom representation
    recover_aaDNA_model(cg_dna, vecXYZ, allBeads)



'''
to do: 
reproduce one of the figures from D. Tree's paper (to verify the model is programmed the way they did)
if the model is fine and good for the length scale of interest then go to the mapping
plotting the energies 9


sassie: mixed monte carlo -> molecular
'''
    # all_atom_pdb = 'trimer_stacked.pdb'

    #s figure()
    #s x = range(iters)
    #s plot(x,rg)
    #s xlabel('iteration #')
    #s ylabel('Rg')

    #test = np.eye(4)
    #s est = np.zeros((10,4))
    #s est[:,1] = range(10)
    #s est[:,3] = 1
    #s est[:,0] = range(10)
    #s est[:,0] *= -1
    #s or i in xrange(4):
    #s        (T, Ti) = move2origin(test[i:])
    #s        test[i:] = np.dot(test[i:],T)
    #s        print 'T:\n', T
    #s        print 'moved2origin:\n', test[i:]
    #s 
    #s        (A, Ai) = align2z(test[i:])
    #s        print 'before align: \n', test[i:]
    #s        test[i:] = np.dot(test[i:],A)
    #s        print 'x A ->\n', A
    #s        print 'aligned: \n', test[i:]
    #s        test[i:] = np.dot(test[i:],Ai)
    #s        print 'un-aligned:  \n', test[i:]
    #s 
    #s        test[i:] = np.dot(test[i:],Ti)
    #s        print 'back2original:\n', test[i:]
    #s 
    #s rint test