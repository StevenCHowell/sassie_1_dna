% $Id: cgDNA_move.py,v 1.3 2013-11-20 16:00:51 schowell Exp $
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

        aa_dna = sasmol.SasMol(0)
        aa_dna.read_pdb(all_atom_pdb)
        natoms = aa_dna.natoms()

#        print 'natoms = ',natoms

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

        coor = cg_dna.coor()[0]

        comment = 'test'
        write_xyz('cg_dna.xyz',coor,comment,frame)
        infile=open('dum.pdb','w')
        for i in xrange(cg_dna.natoms()):
                this_index = index[i]

                sx = cg_dna._coor[frame,i,0]#[:8]
                sy = cg_dna._coor[frame,i,1]#[:8]
                sz = cg_dna._coor[frame,i,2]#[:8]

                infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s\n" % (cg_dna.atom()[i],this_index,cg_dna.name()[i],cg_dna.loc()[i],cg_dna.resname()[i],cg_dna.chain()[i],cg_dna.resid()[i],cg_dna.rescode()[i],sx,sy,sz,cg_dna.occupancy()[i],cg_dna.beta()[i],cg_dna.segname()[i],cg_dna.element()[i],cg_dna.charge()[i]))

        infile.write('END\n')
        infile.close()

        cg_dna.write_pdb("cg_test.pdb",frame,'w')
        cg_dna.write_dcd("cg_test.dcd")

        return cg_dna

def beadRotate(coor3,thetaX,thetaY):

        # add a fourth column of ones to all the coodinates
        (natoms,col) = coor3.shape
        coor4 = np.ones((natoms,4),np.float)
        coor4[:,0:3] = coor3
        # print 'coor4 = ',coor4

        # create the translation-rotation matrix
        # This is intended to be multiplied from the right
        # so as not to require changing the coordinate array.
        tx = thetaX*(np.pi/180.0)
        ty = thetaY*(np.pi/180.0)
        cy = np.cos(ty)
        sy = np.sin(ty)
        cx = np.cos(tx)
        sx = np.sin(tx)
        (x, y, z) = coor3[0,:]
        R = np.eye(4,dtype=np.float)
        (R[0][0], R[0][1], R[0][2]) = (   cy,   0,   -sy)
        (R[1][0], R[1][1], R[1][2]) = (sx*sy,  cx, sx*cy)
        (R[2][0], R[2][1], R[2][2]) = (cx*sy, -sx, cx*cy)
        (R[3][0], R[3][1], R[3][2]) = (-x*cy-y*sx*sy-z*cx*sy +x, -y*cx+z*sx+y, x*sy-y*sx*cy-z*cx*cy+z)
        #print 'R:\n', R

        # this serves to check R
        #T  = np.eye(4,dtype=np.float)
        #Ti = np.eye(4,dtype=np.float) 
        #T[3,0:3] = -coor3[0,:]
        #Ti[3,0:3] = coor3[0,:]
        #Rx = np.eye(4,dtype=np.float)
        #Ry = np.eye(4,dtype=np.float)
        #(Ry[0][0], Ry[0][2]) = (cy, -sy)
        #(Ry[2][0], Ry[2][2]) = (sy,  cy)
        #(Rx[1][1], Rx[1][2]) = ( cx, sx)
        #(Rx[2][1], Rx[2][2]) = (-sx, cx)    
        #Rxy = np.dot(Rx,Ry)
        #print 'Rx:\n', Rx
        #print 'Ry:\n', Ry
        #R_check = np.dot(np.dot(T,Rxy),Ti)
        #print 'R_check:\n', R_check
        
        coor4 = np.dot(coor4,R)

        return coor4[1:,0:3]

def checkU(coor):
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


def dna_mc(nsteps,cg_dna,lp,w):

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

        nbeads = cg_dna.natoms()
        coor = np.copy(cg_dna.coor()[0])

        (u, l) = checkU(coor)

        lpl = lp/l

        # calculate the energy using the current positions
        Ub0 = energyBend(lpl,u)
        Uwca0 = energyWCA(w,coor)
        U_T0 = Ub0 + Uwca0

        for i in xrange(nsteps):

                trial_bead = int((nbeads-1)*random.random())+1
                #print 'trial_bead =', trial_bead

                thetaX = 360*random.random()
                thetaY = 360*random.random()

                # generate a newly rotated model
                #s print 'coor before:\n',coor
                coor[trial_bead:] = beadRotate(coor[trial_bead-1:],thetaX,thetaY)
                #s print 'coor after\n',coor

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
                #if True:
                if test < p:
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

        print '\n\n\n\n'
        
        nsteps = 999
        lp = .0530          # persistence length (actual lp=530A)
        w = 20        # width of dsDNA (actual between 22 & 26) should be < l

        all_atom_pdb = 'dna.pdb'

        cg_dna = make_bead_model(all_atom_pdb)

        print cg_dna.coor

        dna_mc(nsteps,cg_dna,lp,w)

'''
to do: 
reproduce one of the figures from D. Tree's paper (to verify the model is programmed the way they did)
if the model is fine and good for the length scale of interest then go to the mapping
plotting the energies 


sassie: mixed monte carlo -> molecular
'''
