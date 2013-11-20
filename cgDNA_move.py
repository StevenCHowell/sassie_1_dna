import sassie.sasmol.sasmol as sasmol
import numpy,string,os,locale,sys,random,struct

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


def energy():

	return


def rotate():

	return

	
def dna_mc(nsteps,cg_dna):
	
	nbeads = cg_dna.natoms()

	for i in xrange(nsteps):

		trial_bead = int((nbeads-1)*random.random())+1

		coor_terminal = cg_dna.coor()[0][trial_bead]
		coor_base = cg_dna.coor()[0][trial_bead-1]

		uk = coor_terminal-coor_base

#		print uk

	#def rotate(self,frame,axis,theta)
	
	return


def make_bead_model(all_atom_pdb):

	aa_dna = sasmol.SasMol(0)
	aa_dna.read_pdb(all_atom_pdb)
	natoms = aa_dna.natoms()

	#align_dna(aa_dna)

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
	new_coor = numpy.zeros((1,cg_natoms,3),numpy.float)
	print new_coor

	for i in xrange(cg_natoms):
		print "i = ",i
		new_coor[0][i][2] = i*diameter
	
	cg_dna.setCoor(new_coor)
	print new_coor
	print cg_dna.coor()

#	charge = cg_dna.beta()
#	cg_dna.setCharge(charge)

#	print 'loc = ',cg_dna.loc()
#	print 'rescode = ',cg_dna.rescode()
#	print 'charge = ',cg_dna.charge()

#	resid = cg_dna.resid() ; print 'resid = ',resid
#	sresid = numpy.char.mod('%s', resid) ; sresid = sresid.tolist()
#	cg_dna.setResid(sresid)
#	resid = cg_dna.resid() ; print 'resid = ',resid

	print 'type(index) = ',type(cg_dna.index())
	index = cg_dna.index()
#	sindex = str(index)
#	cg_dna.setIndex(sindex)
	
	coor = cg_dna.coor()[0]

	comment = 'test'
	write_xyz('cg_dna.xyz',coor,comment,frame)
	infile=open('dum.pdb','w')
	for i in xrange(cg_dna.natoms()):
		this_index = index[i]

		sx = cg_dna._coor[frame,i,0]#[:8]
		sy = cg_dna._coor[frame,i,1]#[:8]
		sz = cg_dna._coor[frame,i,2]#[:8]

#		print 'test type = ',type(cg_dna._coor[frame][:8])
#		print 'test type = ',type(cg_dna._coor[frame][:8])
#		print 'test type = ',type(cg_dna._coor[frame][:8])
#
#		sx = "{0:8.3f}".format(cg_dna._coor[frame,i,0])[:8]
#		sy = "{0:8.3f}".format(cg_dna._coor[frame,i,1])[:8]
#		sz = "{0:8.3f}".format(cg_dna._coor[frame,i,2])[:8]


		#infile.write("%-6s\n" % (cg_dna.atom()[i]))
		infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s\n" % (cg_dna.atom()[i],this_index,cg_dna.name()[i],cg_dna.loc()[i],cg_dna.resname()[i],cg_dna.chain()[i],cg_dna.resid()[i],cg_dna.rescode()[i],sx,sy,sz,cg_dna.occupancy()[i],cg_dna.beta()[i],cg_dna.segname()[i],cg_dna.element()[i],cg_dna.charge()[i]))
		#infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%2s%2s\n" % (cg_dna._atom[i],this_index,cg_dna._name[i],cg_dna._loc[i],cg_dna._resname[i],cg_dna._chain[i],cg_dna._resid[i],cg_dna._rescode[i],cg_dna._coor[frame,i,0],cg_dna._coor[frame,i,1],cg_dna._coor[frame,i,2],cg_dna._occupancy[i],cg_dna._beta[i],cg_dna._segname[i],cg_dna._element[i],cg_dna._charge[i]))

	infile.write('END\n')
	infile.close()

	cg_dna.write_pdb("cg_test.pdb",frame,'w')
	cg_dna.write_dcd("cg_test.dcd")

	return cg_dna

if __name__ == "__main__":

	nsteps = 10000

	all_atom_pdb = 'dna.pdb'

 	cg_dna = make_bead_model(all_atom_pdb)

	dna_mc(nsteps,cg_dna)


