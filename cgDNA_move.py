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


def energyBend(kT,lp,l,u):
        
        (nu,col) = u.shape
        uk0 = u[:nu-1,:]
        uk1 = u[1:,:]
        res = nu-1-np.sum(uk0*uk1)

        return res*kT*lp/l

def energyWCA(kT,w,coor):

        test = 2.**(1./6.)*w
        (N, col) = coor.shape
        res = 0.0
        for i in xrange(N):
                ri = coor[i,:]
                for j in xrange(i+1,N):
                        rj = coor[j,:]
                        rij = np.sqrt((ri[0]-rj[0])**2. + (ri[1]-rj[1])**2. + (ri[2]-rj[2])**2.)
                        if rij < test:
                                res += (w/rij)**12.-(w/rij)**6.+0.25

        return res*4.*kT

def beadRotate(coor3,thetaX,thetaY):

        # add a fourth column of ones to all the coodinates
        (natoms,col) = coor3.shape
        coor4 = np.ones((natoms,4),np.float)
        coor4[:,0:3] = coor3
        print 'coor4 = ',coor4


        # create the translation matrix
        #s T  = np.eye(4,dtype=np.float) 
        #s Ti = np.eye(4,dtype=np.float) 
        #s T[3,0:3] = -coor3[0,:]
        #s Ti[3,0:3] = coor3[0,:]
        #s print 'T = /n',T
        tx = thetaX*(np.pi/180.0)
        ty = thetaY*(np.pi/180.0)
        cy = np.cos(ty)
        sy = np.sin(ty)
        cx = np.cos(tx)
        sx = np.sin(tx)
        (x, y, z) = coor3[0,:]
        R = np.eye(4,dtype=np.float)
        (R[0][0], R[0][1], R[0][2]) = (   cy,   0,   -sy)
        (R[1][0], R[1][1], R[1][2]) = (sx*sy,  cx, sx*sy)
        (R[2][0], R[2][1], R[2][2]) = (cx*sy, -sx, cx*cy)
        (R[3][0], R[3][1], R[3][2]) = (-x*cy-y*sx*sy-z*cx*sy +x, -y*cx+z*sx+y, x*sy-y*sx*cy-z*cx*cy+z)
        coor4 = np.dot(coor4,R)


        return coor4[:,0:3]

        
def dna_mc(nsteps,cg_dna):

        dcdOutFile = cg_dna.open_dcd_write("cg_dna_moves.dcd")        
        nbeads = cg_dna.natoms()

        for i in xrange(nsteps):

                trial_bead = int((nbeads-1)*random.random())+1

                coor_terminal = cg_dna.coor()[0][trial_bead]
                coor_base = cg_dna.coor()[0][trial_bead-1]

                uk = coor_terminal-coor_base
                #s print uk
                
                thetaX = 360*random.random()
                thetaY = 360*random.random()
                
                cg_dna.coor()[0][trial_bead-1:] = beadRotate(cg_dna.coor()[0][trial_bead-1:],thetaX,thetaY)
                cg_dna.write_dcd_step(dcdOutFile,0,0)

        cg_dna.close_dcd_write(dcdOutFile)

        return


def make_bead_model(all_atom_pdb):

        aa_dna = sasmol.SasMol(0)
        aa_dna.read_pdb(all_atom_pdb)
        natoms = aa_dna.natoms()

#        align_dna(aa_dna)

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
        print 'cg_natoms = ',cg_natoms
        new_coor = np.zeros((1,cg_natoms,3),np.float)
        print new_coor

        for i in xrange(cg_natoms):
                print "i = ",i
                new_coor[0][i][2] = i*diameter
        
        cg_dna.setCoor(new_coor)
        print new_coor
        print cg_dna.coor()

#        charge = cg_dna.beta()
#        cg_dna.setCharge(charge)

#        print 'loc = ',cg_dna.loc()
#        print 'rescode = ',cg_dna.rescode()
#        print 'charge = ',cg_dna.charge()

        print 'type(index) = ',type(cg_dna.index())
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

if __name__ == "__main__":

        nsteps = 10

        all_atom_pdb = 'dna.pdb'

        cg_dna = make_bead_model(all_atom_pdb)

        dna_mc(nsteps,cg_dna)


