import sys
import sassie.sasmol.sasmol as sasmol
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monte_carlo.monomer.dihedral_rotate as dihedral_rotate
import sassie.simulate.monte_carlo.monomer.dihedral_monte_carlo as dihedral_monte_carlo
import sassie.simulate.monte_carlo.monomer.step as step

#       DIHEDRAL_SAMPLING
#
#	09/26/2005 	--	gag-dihedral search		:	jc
#	07/14/2008	--	file management changes		:	jc
#	01/12/2011	--	added sasmol support		:	jc
#	01/15/2011	--	generalized rotation basis	:	jc
#	09/06/2013	--	adapted for affine project	:	jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        DIHEDRAL_SAMPLING is a module that contains the functions
        that are used to generate ensembles of structures by varying
	dihedral angles.

	REFERENCE:

    	J. A. D. MacKerell et al.
    	Journal of Physical Chemistry B,  102  3586-3616  (1998)

    	B. R. Brooks et al.
    	Journal of Computational Chemistry  4  187--217  (1983)

'''
def check_error(error):
	if(len(error)>0):
		print 'error = ',error
		sys.exit()

def initialize_system(m1,molecule_type,basis,numranges,reslow,numcont):

	resid = m1.resid()
	first_last_resid = [resid[0],resid[-1]]

	print 'first_last_resid = ',first_last_resid

	print 'resid = ',resid

	txtOutput=None

	if(molecule_type == 'protein'):
		basis_filter = 'name[i] == "'+basis+'"'
		error,basis_mask = m1.get_subset_mask(basis_filter) ; check_error(error)
		basis_m1=sasmol.SasMol(1)
		error = m1.copy_molecule_using_mask(basis_m1,basis_mask,0) ; check_error(error)
		basis_resname = basis_m1.resname()
		basis_resid = basis_m1.resid()	
		
		print 'reslow = ',reslow
		print 'numcont = ',numcont
		print 'numranges = ',numranges
		print 'basis_resname = ',basis_resname
		print 'basis_resid = ',basis_resid

		respsi=[] ; resphi=[]
		dihedral_energy.protein_initialization(respsi,resphi,basis_resid,basis_resname,numranges,reslow,numcont,first_last_resid,txtOutput)
		dihedral_parameters = [respsi,resphi]
	elif(molecule_type == 'rna'):
		#rna_filter = 'name[i] == "P"'
		#rna_filter = 'name[i] == "O5\'"'
		rna_filter = 'name[i] == "'+basis+'"'
		error,basis_mask = m1.get_subset_mask(rna_filter)
		basis_m1=sasmol.SasMol(1)
		error = m1.copy_molecule_using_mask(basis_m1,basis_mask,0)
		basis_resname = basis_m1.resname()
		basis_resid = basis_m1.resid()	

		print 'sum(basis_mask) rna = ',sum(basis_mask)
		print 'len(resid) rna = ',len(resid)
		resalpha = [] ; resbeta = [] ; resgamma = [] ; resdelta = [] ; resepsilon = [] ; reseta = []
		dihedral_energy.rna_initialization(resalpha,resbeta,resgamma,resdelta,resepsilon,reseta,basis_resid,basis_resname,numranges,reslow,numcont,first_last_resid,txtOutput)
		dihedral_parameters = [resalpha,resbeta,resgamma,resdelta,resepsilon,reseta] 


	print 'dihedral parameters = ',dihedral_parameters


	flexible_residues = dihedral_monte_carlo.get_flexible_residues(numranges,reslow,numcont)

	print 'molecule_type = ', molecule_type

	residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(m1,molecule_type,flexible_residues,txtOutput)

	step_parameters = step.Setup()


	return dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid

if __name__ == '__main__':


	##### REQUIRED INPUT AND SETUP #####

	molecule_type = 'protein'
	temperature = 300.0
	basis = 'CA'
	number_of_ranges = 2 ; numranges=number_of_ranges
	reslow = [2,18]
	numcont = [3,3]
	dtheta = [30.0,30.0]
	nonbondflag = 0
	seed = [1,2]
	
	kb=1.380658E-23 # J/K
	beta=1.0/(temperature*kb)
	
	if(seed[0] == 1):
		from numpy.random import RandomState
		seed_object = RandomState(seed[1])
	else:
		seed_object = -1 


	pdbfile = 'new_test.pdb'
	path = './'

	m1=sasmol.SasMol(0)
	m1.read_pdb(path+pdbfile,saspdbrx_topology=True)
	nf1=m1.number_of_frames()
	
	##### INITIALIZE VARIABLES FOR SEQUENCE THAT WILL BE ROTATED #####


	dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid = initialize_system(m1,molecule_type,basis,numranges,numcont)

	##### LOCAL STUFF AND EXAMPLE USAGE #####

	number_of_trials = 10
	pairdat=['dum',1,1.0]

	coor=m1.coor()

	for i in xrange(number_of_trials):

		print '.', ; sys.stdout.flush()

		vdi,vdf,indices,this_mask=step_parameters.chooser(coor,m1,pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,nonbondflag,first_last_resid,molecule_type,seed_object)

		an=pairdat[0] ; q0=pairdat[1] ; th=pairdat[2]

		print 'angle : q0 : th = ',an,q0,th



