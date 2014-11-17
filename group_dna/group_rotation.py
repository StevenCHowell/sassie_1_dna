import sassie.sasmol.sasmol as sasmol
import sassie.sasmol.sasmath as sasmath
import numpy
import sys,os,string,locale,math
# import sassie.simulate.monte_carlo.ds_dna_monte_carlo.ds_dna_monte_carlo as ddmc
import dna.ds_dna_monte_carlo as ddmc

sys.path.append('./')
import sh_box_of_balls

def make_groups(hybrid,number_of_groups,residues_in_groups):

    frame = 0
    resid = hybrid.resid()
    group_masks = []
    groups = []
    for i in xrange(number_of_groups):
        this_resids = residues_in_groups[i]
        for j in xrange(len(this_resids)):
            if(j==0):
                basis = 'resid[i] == '+str(this_resids[j])+' '			
            else:
                basis += ' or resid[i] == '+str(this_resids[j])

        print '>> creating basis = ',basis
        error,mask = hybrid.get_subset_mask(basis)
        group_masks.append(mask)
        this_group = sasmol.SasMol(0)
        error = hybrid.copy_molecule_using_mask(this_group,mask,frame)
        groups.append(this_group)

    return groups,group_masks


def get_this_residue_local_coordinates_protein(this_group,residue_to_rotate,first_last_resid):

    lcoor = None
    frame = 0

    #### OPEN ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
    #### OPEN ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
    #### OPEN ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS

    if(residue_to_rotate == first_last_resid[0]):
        next_resid = str(residue_to_rotate + 1)       #### OUCH ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
        this_resid = str(residue_to_rotate)

        basis = '(moltype[i] == "protein" and (resid[i] == '+this_resid+' and (name[i] == "N" or name[i] == "CA" or name[i] == "C")) or (resid[i] == '+next_resid+' and name[i] == "N" ) )'	

    elif(residue_to_rotate == first_last_resid[1]):
        previous_resid = str(residue_to_rotate - 1)       #### OUCH ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
        this_resid = str(residue_to_rotate)
        basis = '(moltype[i] == "protein" and (resid[i] == '+this_resid+' and (name[i] == "N" or name[i] == "CA" or name[i] == "C")) or (resid[i] == '+previous_resid+' and name[i] == "C" ) )'	
    else:	
        previous_resid = str(residue_to_rotate - 1)       #### OUCH ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
        next_resid = str(residue_to_rotate + 1)       #### OUCH ... THIS ASSUMES GROUP RESIDS ARE CONTINUOUS
        this_resid = str(residue_to_rotate)
        basis = '(moltype[i] == "protein" and (resid[i] == '+this_resid+' and (name[i] == "N" or name[i] == "CA" or name[i] == "C")) or (resid[i] == '+next_resid+' and name[i] == "N" ) or (resid[i] == '+previous_resid+' and name[i] == "C" ) )'	

    error,mask = this_group.get_subset_mask(basis)

    error,lcoor = this_group.get_coor_using_mask(frame,mask)

    return lcoor[0]

def get_protein_watch_value(this_group,residue_to_rotate,angle,backward):

    if(angle == "phi"):
        if(backward):
            basis = 'name[i] == "CA" and resid[i] == '+str(residue_to_rotate)
        else:
            basis = 'name[i] == "N" and resid[i] == '+str(residue_to_rotate)
    if(angle == "psi"):
        if(backward):
            basis = 'name[i] == "C" and resid[i] == '+str(residue_to_rotate)
        else:
            basis = 'name[i] == "CA" and resid[i] == '+str(residue_to_rotate)

    error,mask = this_group.get_subset_mask(basis)	
    watch_value = (numpy.nonzero(mask*numpy.arange(1,this_group.natoms()+1))[0])[0]

#	print 'resname = ',this_group.resname()[watch_value],'seg = ',this_group.segname()[watch_value],'name = ',this_group.name()[watch_value],'index = ',this_group.index()[watch_value],'resid = ',this_group.resid()[watch_value]

    return watch_value

def rotate_protein_backbone_dihedral(this_group,residue_to_rotate,angle,theta,backward):

#	print '>> rotating protein protein backbone dihedral'
    frame = 0

    coor = this_group.coor() #[frame]
    resid = this_group.resid()

    if not backward:
        first_last_resid = [resid[0],resid[-1]]
    else:
        first_last_resid = [resid[-1],resid[0]]

    lcoor = get_this_residue_local_coordinates_protein(this_group,residue_to_rotate,first_last_resid)

    q0 = residue_to_rotate

    v=numpy.zeros(3,numpy.float) 
    tee=numpy.identity(4,numpy.float)

    if(q0 == first_last_resid[0]):
        n1  = lcoor[0,:]
        ca1 = lcoor[1,:]
        c1  = lcoor[2,:]
        n2  = lcoor[3,:]
    elif(q0 == first_last_resid[1]):
        c0  = lcoor[0,:]
        n1  = lcoor[1,:]
        ca1 = lcoor[2,:]
        c1  = lcoor[3,:]
    else:
        c0  = lcoor[0,:]
        n1  = lcoor[1,:]
        ca1 = lcoor[2,:]
        c1  = lcoor[3,:]
        n2  = lcoor[4,:]

    if(angle == 'phi'):
        v[0]=ca1[0]-n1[0] ; v[1]=ca1[1]-n1[1] ; v[2]=ca1[2]-n1[2]	
        tee[0][3]=n1[0] 
        tee[1][3]=n1[1] 
        tee[2][3]=n1[2]	

        watch_value = get_protein_watch_value(this_group,residue_to_rotate,angle,backward)

    elif(angle == 'psi'):
        v[0]=c1[0]-ca1[0] ; v[1]=c1[1]-ca1[1] ; v[2]=c1[2]-ca1[2]	
        tee[0][3]=ca1[0] 
        tee[1][3]=ca1[1] 
        tee[2][3]=ca1[2]	

        watch_value = get_protein_watch_value(this_group,residue_to_rotate,angle,backward)

    ### fix for first_last case
    #### OPEN
#		if(q0 == first_last_resid[0]):
#			watch_value=indices[2]-1
#		else:
#			watch_value=indices[3]-1

    lv=math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
    nv=v/lv 

    R = box_of_balls.get_rotation_matrix(nv,theta)

    itee = numpy.matrix(tee).I
    ir = R*itee
    br = tee*ir

    natoms=this_group.natoms()

    if(backward):
        qsel = numpy.arange(0,watch_value,1)
    else:
        qsel = numpy.arange(watch_value+1,natoms)

#	this_group.write_pdb('group.pdb',frame,'a')
    newa=[] 

    for i in xrange(natoms):
        if(backward):
            if(i<watch_value):
                ta=coor[frame][i]
                far=ta.tolist()
                far.append(1.0)
                nta=numpy.array(far)
                nta.shape=(-1,4)
                p=br*numpy.matrix(numpy.transpose(nta)) ; p=numpy.array(p)
                newa.append(numpy.transpose(p))
        else:
            if(i>watch_value):
                ta=coor[frame][i]
                far=ta.tolist()
                far.append(1.0)
                nta=numpy.array(far)
                nta.shape=(-1,4)
                p=br*numpy.matrix(numpy.transpose(nta)) ; p=numpy.array(p)
                newa.append(numpy.transpose(p))

    ncoords=numpy.array(newa)
    ncoords.shape=(-1,4)
    temp=numpy.take(ncoords,(0,1,2),-1)

    '''
	if(angle == 'phi'):
		dihedral_angle = sasmath.dihedral_angle(lcoor[0,:],lcoor[1,:],lcoor[2,:],lcoor[3,:])
	elif(angle == 'psi'):
		dihedral_angle = sasmath.dihedral_angle(lcoor[1,:],lcoor[2,:],lcoor[3,:],lcoor[4,:])

	print 'initial_dihedral_angle = ',dihedral_angle
	outfile = open('da.txt','a')
	outfile.write("%lf\n" % (dihedral_angle))
	outfile.close()
	'''

    totsnz1=numpy.zeros(len(temp)*3,numpy.int32)
    lts=0
    for ts in range(len(temp)):
        totsnz1[lts]=(qsel[ts]*3)
        totsnz1[lts+1]=(qsel[ts]*3)+1
        totsnz1[lts+2]=(qsel[ts]*3)+2
        lts=lts+3

    numpy.put(coor[frame],totsnz1,temp)

    #this_group.write_pdb('group.pdb',frame,'a')

    lcoor = get_this_residue_local_coordinates_protein(this_group,residue_to_rotate,first_last_resid)

    '''	
	if(angle == 'phi'):
		dihedral_angle = sasmath.dihedral_angle(lcoor[0,:],lcoor[1,:],lcoor[2,:],lcoor[3,:])
	elif(angle == 'psi'):
		dihedral_angle = sasmath.dihedral_angle(lcoor[1,:],lcoor[2,:],lcoor[3,:],lcoor[4,:])

	print 'final_dihedral_angle = ',dihedral_angle
	'''

    return

def rotate_a_group(this_group,rotate_type,residue_to_rotate,angle,theta,backward):	

    if(rotate_type == 'protein_backbone_dihedral'):
        rotate_protein_backbone_dihedral(this_group,residue_to_rotate,angle,theta,backward)
    if(rotate_type == 'ds_dna'):
        # cg the group to rotate
        # perform the rotation on the residue and following DNA
        # reverse cg the DNA
        
        (d_coor[trial_bead:], xyz[:, trial_bead:], p_coor_rot) = ddmc.beadRotate(
            d_coor[trial_bead-1:], xyz[:, trial_bead-1:], thetaXYZ,
            softrotation, p_coor_rot) 

    return

def hybrid_main(pdb_file_name,number_of_groups,residues_in_groups,rotate_type,group_to_rotate,residue_to_rotate,angle,theta,backward):

    hybrid = sasmol.SasMol(0)
    hybrid.read_pdb(pdb_file_name)

    groups,group_masks = make_groups(hybrid,number_of_groups,residues_in_groups)	

    this_group = groups[group_to_rotate]

    itheta = theta

    for i in xrange(100):
        theta = theta + (5.0 * numpy.pi/180.0)	
        rotate_a_group(this_group,rotate_type,residue_to_rotate,angle,theta,backward)	

    return

def dna_main(pdb_file_name, number_of_groups, residues_in_groups, 
             rotate_type, group_to_rotate, residue_to_rotate, angle, theta,
             backward, dna_segnames, dna_resids, bp_per_bead):



    if rotate_type == 'ds_dna':
        (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, trialbeads, beadgroups, 
         move_masks, all_beads, dna_bead_masks, aa_pgroup_masks, 
         cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask,
         bp_per_bead_array, pkl_file) = ddmc.get_cg_parameters(
             residues_in_groups, dna_resids, dna_segnames, pdb_file_name, 
             pdb_file_name, bp_per_bead)

        # determine which bead has the 'residue_to_rotate' in it
        bead_to_rotate = 9
        
        for i in xrange(100):
            theta = theta + (5.0 * numpy.pi/180.0)    # treating this as thetaX
            thetaXYZ = [theta,0,0]
            
            # rotate that bead
            (cg_dna.coor()[0][bead_to_rotate:], vecXYZ[:,bead_to_rotate:], dummy
             ) = ddmc.beadRotate(cg_dna.coor()[0][bead_to_rotate-1:], 
                                 vecXYZ[:,bead_to_rotate-1:], thetaXYZ,
                                 numpy.zeros((0, 3)) ) 

            #s rotate_dna_group(this_group, rotate_type, residue_to_rotate, angle, theta, backward)	
            
    else:
        mol = sasmol.SasMol(0)
        mol.read_pdb(pdb_file_name)
        groups, group_masks = make_groups(mol, number_of_groups, residues_in_groups)	

        this_group = groups[group_to_rotate]

        itheta = theta

        for i in xrange(100):
            theta = theta + (5.0 * numpy.pi/180.0)	
            rotate_a_group(this_group, rotate_type, residue_to_rotate, angle, theta, backward)	

    return

if __name__ == "__main__":

    # pdb_file_name = 'hybrid.pdb'
    pdb_file_name = 'new_dsDNA60.pdb'
    
    group_to_rotate = 1
    residue_to_rotate = 9
    # angle = 'phi'
    angle = 'psi'
    
    theta = 5.0 ; theta = theta * numpy.pi/180.0
    # backward = False	
    backward = True	

    # rotate_type = 'protein_backbone_dihedral'
    rotate_type = 'ds_dna'
    
    number_of_groups = 5
    # residues_in_groups = [[1,2,3,4,5,6], [7,8,9,10,11,12]]

    # DNA                 chain 1        chain 2 
    # residues_in_groups = [range(1,6)   + range(116,121), 
                          # range(6,11)  + range(111,116), 
                          # range(11,31) + range(91,111),
                          # range(31,56) + range(66,91), 
                          # range(56,61) + range(61,66)] #residues from DNA chain 1
    residues_in_groups = [range(1,6), range(6,11), range(11,31),
                          range(31,56), range(56,61)]
                          
    dna_segnames = ['DNA1','DNA2']
    dna_resids = [[1,60], [120, 61]]
    bp_per_bead = 1
    
    # hybrid_main(pdb_file_name,number_of_groups,residues_in_groups,rotate_type,group_to_rotate,residue_to_rotate,angle,theta,backward)
    dna_main(pdb_file_name, number_of_groups, residues_in_groups, 
             rotate_type, group_to_rotate, residue_to_rotate, angle, theta,
             backward, dna_segnames, dna_resids, bp_per_bead)


'''
ATOM      1 C    MAB     1       0.000   0.000  86.860  0.00  0.00      FAB1 C  
ATOM      2 N    THR     2      -0.346   0.458  61.304  0.00  0.00      LNK1 C  
ATOM      3 CA   THR     2       0.003   0.779  60.103  0.00  0.00      LNK1 C  
ATOM      4 C    THR     2      -0.480   0.006  59.059  0.00  0.00      LNK1 C  
ATOM      5 N    VAL     3      -0.167   0.289  57.944  0.00  0.00      LNK1 C  
ATOM      6 CA   VAL     3      -0.535  -0.333  56.882  0.00  0.00      LNK1 C  
ATOM      7 C    VAL     3       0.032   0.300  55.849  0.00  0.00      LNK1 C  
ATOM      8 N    GLU     4      -0.171  -0.105  54.743  0.00  0.00      LNK1 C  
ATOM      9 CA   GLU     4       0.309   0.406  53.675  0.00  0.00      LNK1 C  
ATOM     10 C    GLU     4      -0.131  -0.310  52.627  0.00  0.00      LNK1 C  
ATOM     11 N    ARG     5       0.215   0.019  51.530  0.00  0.00      LNK1 C  
ATOM     12 CA   ARG     5      -0.124  -0.566  50.437  0.00  0.00      LNK1 C  
ATOM     13 C    ARG     5       0.427   0.041  49.375  0.00  0.00      LNK1 C  
ATOM     14 N    LYS     6       0.195  -0.402  48.286  0.00  0.00      LNK1 C  
ATOM     15 CA   LYS     6       0.645   0.065  47.173  0.00  0.00      LNK1 C  
ATOM     16 C    LYS     6       0.171  -0.691  46.174  0.00  0.00      LNK1 C  
ATOM     17 C    MAB     7      51.666   0.000  69.823  0.00  0.00      FAB2 C  
ATOM     18 N    THR     8      36.187   0.458  49.485  0.00  0.00      LNK2 C  
ATOM     19 CA   THR     8      35.753   0.779  48.312  0.00  0.00      LNK2 C  
ATOM     20 C    THR     8      34.744   0.006  47.760  0.00  0.00      LNK2 C  
ATOM     21 N    VAL     9      34.332   0.289  46.678  0.00  0.00      LNK2 C  
ATOM     22 CA   VAL     9      33.404  -0.333  46.043  0.00  0.00      LNK2 C  
ATOM     23 C    VAL     9      33.246   0.300  44.876  0.00  0.00      LNK2 C  
ATOM     24 N    GLU    10      32.425  -0.105  44.107  0.00  0.00      LNK2 C  
ATOM     25 CA   GLU    10      32.176   0.406  42.963  0.00  0.00      LNK2 C  
ATOM     26 C    GLU    10      31.199  -0.310  42.382  0.00  0.00      LNK2 C  
ATOM     27 N    ARG    11      30.824   0.019  41.295  0.00  0.00      LNK2 C  
ATOM     28 CA   ARG    11      29.902  -0.566  40.618  0.00  0.00      LNK2 C  
ATOM     29 C    ARG    11      29.713   0.041  39.436  0.00  0.00      LNK2 C  
ATOM     30 N    LYS    12      28.878  -0.402  38.699  0.00  0.00      LNK2 C  
ATOM     31 CA   LYS    12      28.578   0.065  37.537  0.00  0.00      LNK2 C  
ATOM     32 C    LYS    12      27.603  -0.691  37.016  0.00  0.00      LNK2 C  
ATOM     33 C    MAB    13       0.000   0.000   0.000  0.00  0.00      FC   C  
END
'''