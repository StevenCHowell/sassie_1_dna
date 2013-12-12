#!/usr/bin/python
# $Id: cgDNA_move.py,v 1.8 2013-12-12 18:25:20 schowell Exp $
import sassie.sasmol.sasmol as sasmol
import numpy as np,string,os,locale,sys,random

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

        print 'natoms = ',natoms

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

def make_model(all_atom_pdb,chain1,chain2,atoms1,atoms2):

        aa_dna = sasmol.SasMol(0)
        aa_dna.read_pdb(all_atom_pdb)
        natoms = aa_dna.natoms()



        print 'natoms = ',natoms

        cg_dna = sasmol.SasMol(0)
        # select the cg beads coordinates from the PDB
        # !!! TODO: need to incorporate user selection of which chains and residues to use as flexible DNA then keep everything eles as rigid bodies !!!
        basis_filter = "(resname[i] == 'ADE' or resname[i] == 'GUA') and (name[i] == 'N1' and (resid[i] < 21) and segname[i]=='DNA1')"
        #s tetramer PDB uses different DNA labelling
        #s basis_filter = "(resname[i] == 'DA' or resname[i] == 'DG')"# and (name[i] == 'N1' and (resid[i] < 21) and segname[i]=='DNA1')"
        error,mask = aa_dna.get_subset_mask(basis_filter) 

        frame = 0
        error = aa_dna.copy_molecule_using_mask(cg_dna,mask,frame) 
        error,coor=aa_dna.get_coor_using_mask(frame,mask)

        cg_dna.setCoor(coor)

        diameter = 6.0*3.4

        cg_natoms = cg_dna.natoms()
        print 'cg_natoms = ',cg_natoms
        new_coor = np.zeros((1,cg_natoms,3),np.float)
#        print new_coor

        for i in xrange(cg_natoms):
#                print "i = ",i
                new_coor[0][i][2] = i*diameter

        cg_dna.setCoor(new_coor)
        print new_coor
        print cg_dna.coor()

#        charge = cg_dna.beta()
#        cg_dna.setCharge(charge)

        #s print 'loc = ',cg_dna.loc()
        #s print 'rescode = ',cg_dna.rescode()
        #s print 'charge = ',cg_dna.charge()

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

        #sh need to generate the local coordinates of all the atoms in each bead relative to the bead's origin
        #sh need to generate the local coordinates of each bead
        vecXYZ = np.zeros((cg_natoms*3,3))
        vecXYZ[0:cg_natoms] = [1,0,0]
        vecXYZ[cg_natoms:2*cg_natoms] = [0,1,0]
        vecXYZ[2*cg_natoms:3*cg_natoms] = [0,0,1]
        return (cg_dna, vecXYZ)

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
        Axz  = np.eye(4,dtype=np.float)
        Az   = np.eye(4,dtype=np.float)
        assert all(coor4[0] == [0., 0., 0., 1.,]), "coordinates passed to align2z were not translated to the origin"

        if coor4.shape > (1, 4):
                small = 1E-14 # small limit to make sure not dividing by zero
                
                # (u, v, w) pointing from bead0 to bead1
                (u, v, w) = coor4[1,0:3]
                print '(u, v, w) = ', (u,v,w)
                
                # align to the x-z plane
                d1 = np.sqrt(u**2+v**2)
                if v > small:
                        # print '(u, v, d1) =', (u, v, d1)
                        (Axz[0][0], Axz[0][1]) = (u/d1, -v/d1)
                        (Axz[1][0], Axz[1][1]) = (v/d1,  u/d1)  
                        print 'Axz= \n', Axz
                else:
                        print 'already aligned to xz-plane'

                        # align to the z-axis
                d2 = np.sqrt(u**2+v**2+w**2)
                if d1 > small:
                        (Az[0][0], Az[0][2]) = (w/d2,  d1/d2)
                        (Az[2][0], Az[2][2]) = (-d1/d2, w/d2)
                        print 'Az= \n', Az
                else:
                        print 'already aligned to z-axis'                        
        else:
                print 'no point to align'

        A = np.dot(Axz,Az)
        Ai = np.dot(Az.transpose(),Axz.transpose())
        # print 'coor passed to A, Ai:\n',coor4
        #s print '[coor] x [Axz]:\n',np.dot(coor4,Axz)
        # print 'AxzI test:\n',np.dot(np.dot(coor4,Axz),Axz.transpose())
        #s print "Axz = \n", Axz
        #s print "Az = \n", Az
        #s print "these 3 should be identities:"
        #s print np.dot(Az,Az.transpose())
        #s print np.dot(Axz,Axz.transpose())
        #s print np.dot(A,Ai)
        return (A, Ai)
        

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

        coor4 = np.ones((natoms,4),np.float)
        X = np.copy(coor4)
        Y = np.copy(coor4)
        Z = np.copy(coor4)
        coor4[:,0:3] = coor3
        X[:,0:3] = vecX
        Y[:,0:3] = vecY
        Z[:,0:3] = vecZ
        # print 'coor4 = ',coor4

        # create the translation-rotation matrix
        # This is intended to be multiplied from the right (unlike standard matrix multiplication)
        # so as not to require transing the coordinate vectors.
        tx = thetas[0]*(np.pi/180.0)
        ty = thetas[1]*(np.pi/180.0)
        tz = thetas[2]*(np.pi/180.0)
        cx = np.cos(tx)
        sx = np.sin(tx)
        cy = np.cos(ty)
        sy = np.sin(ty)
        cz = np.cos(tz)
        sz = np.sin(tz)
        #(x, y, z) = coor3[0,:]
        #print 'HERE ---->',  (x, y, z)
        
        # initialize the transformation pieces
        Rx   = np.eye(4,dtype=np.float)
        Ry   = np.eye(4,dtype=np.float)
        Rz   = np.eye(4,dtype=np.float)

        (Rx[1][1], Rx[1][2]) = ( cx, sx)
        (Rx[2][1], Rx[2][2]) = (-sx, cx)  

        (Ry[0][0], Ry[0][2]) = (cy, -sy)
        (Ry[2][0], Ry[2][2]) = (sy,  cy)

        (Rz[0][0], Rz[0][1]) = ( cz, sz)
        (Rz[1][0], Rz[1][1]) = (-sz, cz)  
        # print 'Rx:\n', Rx
        # print 'Ry:\n', Ry
        # print 'Rz:\n', Rz

        print 'original coor:\n', coor4
        (T0, Ti0) = move2origin(coor4)
        coor4 = np.dot(coor4,T0)  # move to origin
        print 'moved2origin coor:\n', coor4
                
        (A, Ai) = align2z(coor4)
        coor4 = np.dot(coor4,A) # align to z-axis 
        print 'aligned coor:\n', coor4

        Rxyz = np.dot(np.dot(Rx, Ry), Rz)
        #s print "Rxyz = \n", Rxyz
        coor4 = np.dot(coor4,Rxyz) # rotate about first angle
        print 'step 0 rotated coor:\n', coor4

        # repeat rotation for softening the bend
        for i in xrange(1,nSoft):
                (T, Ti) = move2origin(coor4[i:])
                coor4[i:] = np.dot(coor4[i:],T)  # move to origin
                print 'moved2origin coor:\n',coor4

                coor4[i:] = np.dot(coor4[i:],Rxyz)
                print 'step %d' %i,'rotated coor:\n',coor4

                coor4[i:] = np.dot(coor4[i:],Ti) # return to original position
                print 'returned from origin coor:\n',coor4
                # the coarse grained beads local coordinates should not be translated, only rotated
                X[i:] = np.dot(X[i:],Rxyz)
                Y[i:] = np.dot(Y[i:],Rxyz)
                Z[i:] = np.dot(Z[i:],Rxyz)

        coor4 = np.dot(coor4,Ai)
        print 'un-aligned:\n',coor4
        coor4 = np.dot(coor4,Ti0)
        print 'returned from origin coor:\n',coor4
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
        l = np.mean(lu[:])                          # average distance btwn bead 

        #print '\n\ncoor:\n', coor
        #print 'u-vectors:\n', u
        #print 'magnitude of u-vectors:\n', lu
        #print 'average magnitude of u =', l
        #print '\n\n'

        test = np.abs(lu-l) > erLim  # check if any of the l-lengths are different
        if test.any():
                print 'ERROR: the beads are not uniformly spaced'
                print 'u = \n',u
                print test
                # print 'difference = ', u[:,2] - l
                print 'distance = ', lu
        
        return (u, l)

def energyBend(lpl,u):
        ''' 
        this function finds the E/kT for N beads connected by N-1 rods
        lpl is the persistence length divided by the rod length: lp/l
        u contains the N-1 vectors representing each rod
        '''

        (nu,col) = u.shape          # nu = nbeads - 1
        uk0 = u[:nu-1,:]            # define u_k
        uk1 = u[1:,:]               # define u_k+1
        res = nu-1-np.sum(uk0*uk1)  # calculate the sum of their products

        return res*lpl

def energyWCA(w,coor):
        '''
        this function finds the Weeks-Chandler-Anderson repulsion E/kt for N
        beads connected by N-1 rods together representing a worm-like chain.
        w is the effective width of the chain
        coor contains the xyz-coordinates of the beads
        '''
        
        test = 2.**(1./6.)*w
        (N, col) = coor.shape
        res = 0.0
        for i in xrange(N):
                ri = coor[i,:]
                for j in xrange(i+1,N):
                        rj = coor[j,:]
                        rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)
                        # print '(2^(1/6)*w, rij) = ', (test, rij)
                        if rij < test:
                                res += (w/rij)**12.-(w/rij)**6.+0.25

        #s print 'U_wca =', res*4
        return res*4


def dna_mc(nsteps,cg_dna,vecXYZ,lp,w):
        '''
        this function perform nsteps Monte-Carlo moves on the cg_dna
        '''
        nSoft = 3  # number of beads to seperate the move over

        dcdOutFile = cg_dna.open_dcd_write("cg_dna_moves.dcd")        
        import time
        timestr = time.strftime("%y%m%d")
        fName = timestr + 'dnaMoves.o'
        print '\n' + fName + ' contains result paramaters'
        #s if os.path.exists(fName):
        #s         header = ''
        #s else:
        # (dU, Ub1-Ub0, Uwca1-Uwca0, p, test, w, lp, trial_bead, thetaX, thetaY)
        header = '# dU\t\t dU_b\t\t dU_WCA\t\t prob\t random\t w\t lp\t bead\t thetaX\t thetaY\n'

        outData = open(fName,'a')
        outData.write(header)

        nbeads = cg_dna.natoms()  # need to change this so there could be components that are not the beads
        coor = np.copy(cg_dna.coor()[0])

        (u, l) = checkU(coor) # get the vectors pointing from bead to bead and make sure they are equidistant

        lpl = lp/l  # setup the presistence length paramater

        # calculate the energy of the starting positions
        Ub0 = energyBend(lpl,u)
        Uwca0 = energyWCA(w,coor)
        U_T0 = Ub0 + Uwca0

        # split vecXYZ into three matrices
        cg_natoms = cg_dna.natoms()
        vecX = vecXYZ[0:cg_natoms]
        vecY = vecXYZ[cg_natoms:2*cg_natoms]
        vecZ = vecXYZ[2*cg_natoms:3*cg_natoms]

        for i in xrange(nsteps):

                trial_bead = int((nbeads-1)*random.random())+1
                print 'trial_bead =', trial_bead

                thetaX = 180*random.random() - 90
                thetaY = (180*random.random() - 90) * 0
                thetaZ = (15*random.random() - 7.5) * 0# want to make this about +/-5 degrees
                thetaXYZ = [thetaX/nSoft, thetaY/nSoft, thetaZ/nSoft]
                
                # generate a newly rotated model
                # print vecX[trial_bead-1:]
                # print vecY[trial_bead-1:]
                # print vecZ[trial_bead-1:]

                #print 'coor before:\n',coor
                
                #s tried to do x, y, then z but had problems with it
                #s for i in xrange(3):
                #s         thetas = np.zeros(3)
                #s         thetas[i] = thetaXYZ[i]
                #s         print 'new thetas:', thetas

                (coor[trial_bead:],vecX[trial_bead:],vecY[trial_bead:],vecZ[trial_bead:]) = beadRotate(coor[trial_bead-1:],vecX[trial_bead-1:],vecY[trial_bead-1:],vecZ[trial_bead-1:],thetaXYZ,nSoft)
                print 'coor after\n',coor

                # calculate the change in energy (dU) and boltzman factor (p) for the new model
                (u, l) = checkU(coor)
                Ub1 = energyBend(lpl,u)
                Uwca1 = energyWCA(w,coor)
                U_T1 =  Ub1 + Uwca1
                dU = U_T1 - U_T0
                p = np.exp(-dU)
                test = random.random()

                # output the results to a text file
                outData.write("%1.3e\t %1.3e\t %1.3e\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %d\t %0.1f\t %0.1f\n" %(dU, Ub1-Ub0, Uwca1-Uwca0, p, test, w, lp, trial_bead, thetaX, thetaY) )
                # print '(dU,dUb,dUwca,p,rand) =', (dU,Ub1-Ub0, Uwca1-Uwca0, p, test)

                # if accepted write new coordinates, else write old again
                if True:
                #if test < p:
                        print 'wrote new dcd frame (end of loop',i,' trial_bead=',trial_bead,')'+' accepted new configuration'
                        cg_dna.coor()[0] = np.copy(coor)
                else :
                        print 'wrote new dcd frame (end of loop',i,' trial_bead=',trial_bead,')'
                        coor = np.copy(cg_dna.coor()[0])   # reset the coordinates

                cg_dna.write_dcd_step(dcdOutFile,0,0)
                
        # close output files
        cg_dna.close_dcd_write(dcdOutFile)
        outData.close()

        return


if __name__ == "__main__":

        print '\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'

        nsteps = 10
        lp = .0530          # persistence length (actual lp=530A)
        w = 20        # width of dsDNA (actual between 22 & 26) should be < l

        all_atom_pdb = 'dna.pdb'
        #s all_atom_pdb = '1zbb_tetra.pdb'


        (cg_dna, vecXYZ) = make_model(all_atom_pdb, 1, 2, range(12), range(12))
        print cg_dna.coor
        dna_mc(nsteps,cg_dna,vecXYZ,lp,w)
        
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

'''
to do: 
reproduce one of the figures from D. Tree's paper (to verify the model is programmed the way they did)
if the model is fine and good for the length scale of interest then go to the mapping
plotting the energies 


sassie: mixed monte carlo -> molecular
'''
