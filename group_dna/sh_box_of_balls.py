import sassie.sasmol.sasmol as sasmol
import numpy,math,random
import os,sys,locale
import sassie.simulate.rigid_body.two_body_grid.poverlap as poverlap 
import sassie.sasmol.sasmath as sasmath
import sassie.simulate.monte_carlo.monomer.dihedral_rotate as dihedral_rotate

import pdb

# sys.path.append('/Users/curtisj/Desktop/august_sassie_development/svn_utk/sassie_1.0/new_modules/mmc')
# sys.path.append('/Users/curtisj/sassie_full_svn/sassie_1.0/new_modules/mmc')
sys.path.append('/home/myPrograms/svn_utk/new_modules/mmc')

import setup_mmc
import move
sys.path.append('./')
import initialize_rotation 
import group_rotation

import test_rotation_matrix

def check_error(error):
    if(len(error)>0):
        print error
        sys.exit()

def check_ball_overlap(single_molecule,this_group,group_to_rotate,cutoff):


    maxdist = 86.859784 + 3.0	
    check = 1

    frame = 0

    #### HARDWIRED

    fc_basis = 'segname[i] == "FC"'
    fab1_basis = 'segname[i] == "FAB1"'
    fab2_basis = 'segname[i] == "FAB2"'

    error,fc_mask = single_molecule.get_subset_mask(fc_basis) ; check_error(error)
    error,fab1_mask = single_molecule.get_subset_mask(fab1_basis) ; check_error(error)
    error,fab2_mask = single_molecule.get_subset_mask(fab2_basis) ; check_error(error)

    error,fc_coor = single_molecule.get_coor_using_mask(frame,fc_mask) ; check_error(error)

    #fc_coor = single_molecule.coor()[0][32][:]
    if(group_to_rotate == 0):
        fab1_coor = this_group.coor()[0][0][:]
        fab2_coor = single_molecule.coor()[0][16][:]
    else:	
        fab1_coor = single_molecule.coor()[0][0][:]
        fab2_coor = this_group.coor()[0][0][:]

#	print '>>> DEBUG CHECK \n\nfab1_coor = ',fab1_coor,'\n\n'
#	print '>>> DEBUG CHECK \n\nfab2_coor = ',fab2_coor,'\n\n'
#	print '>>> DEBUG CHECK \n\nfc_coor = ',fc_coor,'\n\n'

    dist1 = math.sqrt(numpy.sum((fc_coor-fab1_coor)**2.0))
    dist2 = math.sqrt(numpy.sum((fc_coor-fab2_coor)**2.0))
    dist3 = math.sqrt(numpy.sum((fab1_coor-fab2_coor)**2.0))


    if(dist1 > cutoff and dist2 > cutoff and dist3 > cutoff):
        check=0


    #### OPEN
    #### OPEN
    #### OPEN

    ### THIS SHOULD NOT BE NEEDED --->>>> THERE IS A BUG IN ROTATION!!!!
    ### THIS SHOULD NOT BE NEEDED --->>>> THERE IS A BUG IN ROTATION!!!!
    ### THIS SHOULD NOT BE NEEDED --->>>> THERE IS A BUG IN ROTATION!!!!

    if((check == 0) and (dist1 > maxdist or dist2 > maxdist)):
        check=1	

    return check

def rotate_single_mab(single_molecule,trial_molecule,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff):

    frame = 0

    seed = [1,2]   # if you want a set SEED
    seed = [0,8]  

    if(seed[0] == 1):
        from numpy.random import RandomState
        seed_object = RandomState(seed[1])
    else:
        seed_object = -1
    nonbondflag = False

    backward = True
    rotate_type = 'protein_backbone_dihedral'

    molecule_type = rotation_variables[0]
    temperature = rotation_variables[1]
    basis = rotation_variables[2]
    number_of_ranges = rotation_variables[3]
    reslow = rotation_variables[4]
    numcont = rotation_variables[5]
    dtheta = rotation_variables[6]

    kb=1.380658E-23 # J/K
    beta=1.0/(temperature*kb)

    pairdat=['dum',1,1.0]

    check = 1
    all_atoms_basis = 'resid[i] < 100'
    error,all_atoms_mask = single_molecule.get_subset_mask(all_atoms_basis); check_error(error)

    while(check == 1):

        print ':', ; sys.stdout.flush()
        coor=single_molecule.coor()

        vdi,vdf,indices,this_mask=step_parameters.chooser(coor,trial_molecule,pairdat,dtheta,number_of_ranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,nonbondflag,first_last_resid,molecule_type,seed_object)

        angle = pairdat[0] ; residue_to_rotate = pairdat[1] ; theta = pairdat[2]

        #### HARDWIRED
        if(residue_to_rotate < 7):
            this_group = groups[0]
            this_mask =  group_masks[0]
            group_to_rotate = 0
        else:
            this_group = groups[1]
            this_mask =  group_masks[1]
            group_to_rotate = 1

        group_rotation.rotate_a_group(this_group,rotate_type,residue_to_rotate,angle,theta,backward)
        check = check_ball_overlap(trial_molecule,this_group,group_to_rotate,cutoff)

        if(check == 1):	
            trial_molecule.set_coor_using_mask(single_molecule,frame,all_atoms_mask)	
            error = single_molecule.copy_molecule_using_mask(this_group,this_mask,frame)

    if(check == 0):
        trial_molecule.set_coor_using_mask(this_group,frame,this_mask)	
        error = single_molecule.copy_molecule_using_mask(this_group,this_mask,frame)

    return



def simulate_single_mab(single_molecule,number_of_steps,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff):

    frame = 0

    seed = [1,2]
    seed = [0,2]

    if(seed[0] == 1):
        from numpy.random import RandomState
        seed_object = RandomState(seed[1])
    else:
        seed_object = -1
    nonbondflag = False

    backward = True
    rotate_type = 'protein_backbone_dihedral'

    molecule_type = rotation_variables[0]
    temperature = rotation_variables[1]
    basis = rotation_variables[2]
    number_of_ranges = rotation_variables[3]
    reslow = rotation_variables[4]
    numcont = rotation_variables[5]
    dtheta = rotation_variables[6]

    kb=1.380658E-23 # J/K
    beta=1.0/(temperature*kb)

    pairdat=['dum',1,1.0]

    dcdoutfile=single_molecule.open_dcd_write('adum.dcd')
    coor=single_molecule.coor()
    taccepted=0 ; nsteps=0 ; count = 0

    ### debugging

    da_outfile = open('dangles.txt','w')

    ### end debugging

    for i in xrange(number_of_steps):

        print '.', ; sys.stdout.flush()

        vdi,vdf,indices,this_mask=step_parameters.chooser(coor,single_molecule,pairdat,dtheta,number_of_ranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,nonbondflag,first_last_resid,molecule_type,seed_object)

        angle = pairdat[0] ; residue_to_rotate = pairdat[1] ; theta = pairdat[2]

        ### debugging

        #angle = 'psi'
        #residue_to_rotate = 5
        #theta = 5.0 * numpy.pi/180.0

        ### end debugging

        #### HARDWIRED
        if(residue_to_rotate < 7):
            this_group = groups[0]
            this_mask =  group_masks[0]
            group_to_rotate = 0
        else:
            this_group = groups[1]
            this_mask =  group_masks[1]
            group_to_rotate = 1

        group_rotation.rotate_a_group(this_group,rotate_type,residue_to_rotate,angle,theta,backward)
        check = check_ball_overlap(single_molecule,this_group,group_to_rotate,cutoff)

        if(check == 0):

            ### debugging

            if(angle == 'psi'):

                basis = '(resid[i] == '+str(residue_to_rotate)+' and (name[i] == "N" or name[i] == "CA" or name[i] == "C")) or (resid[i] == '+str(residue_to_rotate+1)+' and (name[i] == "C") )'
                error,mask = single_molecule.get_subset_mask(basis) ; check_error(error)


                error,temp_coor = single_molecule.get_coor_using_mask(frame,mask) ; check_error(error)
                lcoor = temp_coor[0]	
                d_angle = sasmath.dihedral_angle(lcoor[0],lcoor[1],lcoor[2],lcoor[3])
                print 'initial dihedral angle = ',d_angle
                da_outfile.write("%lf\n" % (d_angle)) ; da_outfile.flush()

            ### end debugging

            single_molecule.set_coor_using_mask(this_group,frame,this_mask)	

            ### debugging

            if(angle == 'psi'):
                error,temp_coor = single_molecule.get_coor_using_mask(frame,mask) ; check_error(error)
                lcoor = temp_coor[0]	
                d_angle = sasmath.dihedral_angle(lcoor[0],lcoor[1],lcoor[2],lcoor[3])
                print 'final dihedral angle = ',d_angle

                port = 54321
                flag = 0
                single_molecule.send_coordinates_to_vmd(port,flag)	

            ### end debugging

            single_molecule.write_dcd_step(dcdoutfile,0,taccepted+1)
            taccepted += 1
        else:
#			if(count == 10):
#				print '.',
#				dcdfile2 = dum_molecule.open_dcd_read('adum.dcd')
#				nf=dcdfile2[2]
#				my_frame = random.randint(0,nf)	
#				for j in xrange(my_frame):
#					dum_molecule.read_dcd_step(dcdfile2,my_frame)
#				dum_coor = dum_molecule.coor()[0]
#				single_molecule.coor()[0] = dum_coor
#				count = 0
            count += 1

        nsteps += 1


    return


def get_coords_from_argon(argon_file_name):

    argon = sasmol.SasMol(0)
    argon.read_pdb(argon_file_name)
    com_coor = argon.coor()[0]

    return com_coor

def get_com_coords_and_build_box(single_molecule,coor_start_type,number_of_molecules,argon_file_name,argon_box_length,length_scale_factor,mab_diameter):

    single_molecule.center(0)
    mm1=single_molecule.calcminmax()

    dx=math.fabs(mm1[0][0]-mm1[1][0])
    dy=math.fabs(mm1[0][1]-mm1[1][1])
    dz=math.fabs(mm1[0][2]-mm1[1][2])

    max_diameter_molecule = max(dx,dy,dz) ; print 'maximum molecular diameter = ',max_diameter_molecule

    if(mab_diameter > 0.0):	
        min_diameter_molecule = mab_diameter
    else:
        min_diameter_molecule = min(dx,dy,dz) ; print 'minimum molecular diameter = ',min_diameter_molecule

    print 'max diameter = ',max_diameter_molecule
    print 'min diameter = ',min_diameter_molecule

    if(coor_start_type==1):
        com_coor,boxl = setup_mmc.makefcc(2.0*max_diameter_molecule,number_of_molecules)

    elif(coor_start_type==2):
        com_coor = get_coords_from_argon(argon_file_name)
        boxl = argon_box_length*length_scale_factor	
        com_coor = com_coor*length_scale_factor
    else:
        print 'incorrect start up type indicated ',start_type
        print '# 1--> make fcc coords, 2--> restart from PDB coords'
        print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

    return com_coor,boxl,max_diameter_molecule

def get_segname(i,molecule):

    if(i<10):
        this_segment = '000'+str(i)     
    elif(i<100):
        this_segment = '00'+str(i)     
    elif(i<1000):
        this_segment = '0'+str(i)       
    else:
        this_segment = ''+str(i)       

    segname = []

    for j in xrange(molecule.natoms()):
        segname.append(this_segment)

    molecule.setSegname(segname)

    return 


def duplicate_molecule_and_assign_coor(all_molecules,single_molecule,com_coor,initial_pdb,number_of_steps,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff):

    frame = 0

    number_of_duplicates = len(com_coor)
    print ' >> making ',number_of_duplicates,' molecular copies'

    bad_overlap = 0

    ball_basis = 'segname[i] == "FAB1" or segname[i] == "FAB2" or segname[i] == "FC"'
    ball_basis = 'resid[i] == 1 or resid[i] == 7 or resid[i] == 13'
    error,ball_mask = single_molecule.get_subset_mask(ball_basis)

    original_single_molecule = sasmol.SasMol(0)
    all_atoms_basis = 'resid[i] < 100'
    error,all_atoms_mask = single_molecule.get_subset_mask(all_atoms_basis); check_error(error)

    error = single_molecule.copy_molecule_using_mask(original_single_molecule,all_atoms_mask,frame) ; check_error(error) 	

    growing_molecule = sasmol.SasMol(0)

    ntr = 0
    max_count = 2E2

    number_accepted_structures = 0

    dcdoutfile2=single_molecule.open_dcd_write('sdum.dcd')

    for i in xrange(number_of_duplicates):
    #for i in xrange(400):
        print i, ; sys.stdout.flush()	

        ### create a new trial molecule object

        trial_molecule = sasmol.SasMol(0)

        ### get the COM position that we want to fill

        this_com = [com_coor[i]]
        ensemble_count = 0

        check = 1
        if(i>0):

            ### note that the check_error(error) merely prints out the error and quits program

            ### grab the current mask for the balls-only from the complete system

            error,growing_mask = growing_molecule.get_subset_mask(ball_basis) ; check_error(error)

            ### get the coordinates of the balls-only from the complete system

            error,coor1 = growing_molecule.get_coor_using_mask(frame,growing_mask) ; check_error(error)

            count = 0  ; rotate_count = 0

            ### look for a new single molecule to add at this_com with a random rotation (x,y,z) about com
            ###		and a random internal rotations of Fab arms

            ### check = 1 means there is overlap

            while(check==1):
                print '.', ; sys.stdout.flush()
                save_coor = True
                ### assign coordinates from original single molecule to trial_molecule where the coordinates
                ###	are randomly rotated about an arbitrary axis (x,y,or z) and the COM is moved to this_com


                #### SHOULD MERELY COPY SINGLE MOLECULE TO TRIAL MOLECULE ONCE THEN REASSIGN COORDS AND ROTATE HERE

                zero = numpy.array([0.0,0.0,0.0])	
                #error = single_molecule.duplicate_molecule(trial_molecule,1,frame,zero) ; check_error(error)
                error = single_molecule.copy_molecule_using_mask(trial_molecule,all_atoms_mask,frame) ; check_error(error)

                rotate_single_mab(single_molecule,trial_molecule,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff)

                axis_int = random.randint(1,3)
                if(axis_int == 1):
                    axis = 'x'
                elif(axis_int == 2):
                    axis = 'y'
                else:
                    axis = 'z'

                trial_molecule.write_dcd_step(dcdoutfile2,0,number_accepted_structures+1)

                rigid_theta = 360.0*random.random()

                trial_molecule.rotate(frame,axis,rigid_theta)

                trial_molecule.moveto(frame,this_com[0])

                ### get the coordinates of the balls-only from the trial_molecule

                error,coor2 = trial_molecule.get_coor_using_mask(frame,ball_mask) ; check_error(error)

                ### check for overlap of the balls between the complete system and the trial_molecule

                check = poverlap.faboverlap(coor1[0],coor2[0],cutoff)
                count+=1


                if(check == 1 and count == max_count):
                    print '>> need to translate and rotate a ball' ; ntr+=1
                    save_coor = False
                    check=0
                    count = 0
                elif(check == 0):
                    number_accepted_structures += 1

                ### the trial_molecule at COM = this_com does not overlap with complete system
                ### 	will break out of while loop 
                ###

                ###


            ### you have left the while loop with a new molecule that does not overlap with the system
            ### 	set the segname and merge the trial_molecule with the last_molecule to form the 
            ###	new growing_molecule
            if save_coor:
                get_segname(number_accepted_structures,trial_molecule)
                growing_molecule.merge_two_molecules(last_molecule,trial_molecule)	

        else:
            error = single_molecule.duplicate_molecule(trial_molecule,1,frame,this_com)
            check_error(error)
            growing_molecule = trial_molecule
            dcdoutfile = growing_molecule.open_dcd_write('adum.dcd')

        print 'bad_overlap = ',bad_overlap,'\tcheck = ',check

        ### set the last_molecule to the current growing_molecule

        last_molecule = growing_molecule	

    #	growing_molecule.write_dcd_step(dcdoutfile,0,i+1)

    growing_molecule.write_pdb("adum.pdb",frame,"w")

    print 'bad overlap = ',bad_overlap
    print 'need to rotate = ',ntr

    number_atoms = all_molecules.natoms()

    print ' >> all natoms = ',number_atoms

    return

def get_pdb_values(m1,natoms):

    atom=[] ; index=[] ; name=[] ; loc=[] ; resname=[] ; chain=[] ; resid=[] ; rescode=[]
    x=[] ; y=[] ; z=[]
    occupancy=[] ; beta=[] ; segname=[] ; element=[] ; charge=[] ; moltype=[]	

    for i in xrange(natoms):
        atom.append("ATOM  ")
        index.append(i+1)
        name.append("C")
        loc.append(" ")
        resname.append("MAB")
        chain.append(" ")
        resid.append(i+1)
        rescode.append(" ")
        occupancy.append("  0.00")	
        beta.append("  0.00")	
        segname.append("MAB1")
        element.append("C")
        charge.append("  ")
        moltype.append("other")
    m1.setAtom(atom) ; m1.setIndex(index) ; m1.setName(name) ; m1.setLoc(loc) ; m1.setResname(resname)
    m1.setChain(chain) ; m1.setResid(resid) ; m1.setRescode(rescode) ; m1.setOccupancy(occupancy)
    m1.setBeta(beta) ; m1.setSegname(segname) ; m1.setElement(element) ; m1.setCharge(charge)
    m1.setMoltype(moltype)
    m1.setNatoms(natoms)

    return

def build_ball_coordinates(ball_diameter,fc_to_fab_vector_length,miniumum_base_angle):

    coor = numpy.zeros((1,3,3),numpy.float)
    coor_1 = numpy.zeros((1,3),numpy.float)
    coor_2 = numpy.zeros((1,3),numpy.float)
    coor_3 = numpy.zeros((1,3),numpy.float)

    coor_2[0][2] = fc_to_fab_vector_length

    m1 = sasmol.SasMol(0)
    dum = numpy.copy(coor_2)
    m1.setCoor(dum)
    frame=0
    m1.rotate(frame,'y',minimum_base_angle)

    coor_3 = m1.coor()

    print 'coor_1 = ',coor_1[0]
    print 'coor_2 = ',coor_2[0]
    print 'coor_3 = ',coor_3[0]

    dist_1_2 = numpy.sqrt(numpy.sum((coor_1[0]-coor_2[0])**2.0))
    dist_1_3 = numpy.sqrt(numpy.sum((coor_1[0]-coor_3[0])**2.0))
    dist_2_3 = numpy.sqrt(numpy.sum((coor_2[0]-coor_3[0])**2.0))

    print 'dist_1_2 = ',dist_1_2
    print 'dist_1_3 = ',dist_1_3
    print 'dist_2_3 = ',dist_2_3
    print 'ball_diameter = ',ball_diameter

    coor[0][0] = coor_1[0]
    coor[0][1] = coor_2[0]
    coor[0][2] = coor_3[0]

    get_pdb_values(m1,3)
    m1.setCoor(coor)
    m1.write_pdb('three_ball.pdb',frame,'w')

    return m1

def get_rotation_matrix(input_vector,theta):
    '''
    returns rotation vector to rotate about the input_vector by theta (radians)
    '''
    c = numpy.cos(theta) ; s = numpy.sin(theta) ; t = 1.0 - c
    x = input_vector[0] ; y = input_vector[1] ; z = input_vector[2]
    R = numpy.zeros((4,4),numpy.float)
    R[0][0] = t*x*x + c
    R[0][1] = t*x*y - s*z
    R[0][2] = t*x*z + s*y
    R[1][0] = t*x*y + s*z
    R[1][1] = t*y*y + c
    R[1][2] = t*y*z - s*x
    R[2][0] = t*x*z - s*y
    R[2][1] = t*y*z + s*x
    R[2][2] = t*z*z + c
    R[3][3] = 1.0

    R=numpy.matrix(R)

    test_rotation_matrix.evaluate_rotation_matrix(R,theta)

    return R

def get_alignment_rotation_matrix(b,a):
    '''
    returns general rotation matrix to rotate vector b onto vector a
    '''
    mag_a = numpy.sqrt(numpy.sum(a*a)) 
    mag_b = numpy.sqrt(numpy.sum(b*b)) 
    a = a/mag_a
    b = b/mag_b
    a_dot_b = numpy.dot(a,b)
    theta = numpy.arccos(a_dot_b)	
    b_cross_a = sasmath.cross_product(b,a)
    R = get_rotation_matrix(b_cross_a,theta)

    return R

def rotate_coordinates(b,R,linker):

    '''
	for an end-point vector b (origin at zero) and a pre-determined rotation matrix
	this routine rotates the coordinates from frame=0 of object linker using R & b
	'''

    frame = 0
    v=numpy.zeros(3,numpy.float) 
    tee=numpy.identity(4,numpy.float)
    tee[0][3] = b[0]
    tee[1][3] = b[1]
    tee[2][3] = b[2]

    itee = numpy.matrix(tee).I
    ir = R*itee
    br = tee*ir

    natoms=linker.natoms()
    linker_coor = linker.coor()
    newa=[]
    #for i in xrange(natoms-1):
    for i in xrange(natoms):
        ta=linker_coor[frame][i]
        far=ta.tolist()
        far.append(1.0)
        nta=numpy.array(far)
        nta.shape=(-1,4)
        p=br*numpy.matrix(numpy.transpose(nta)) ; p=numpy.array(p)
        newa.append(numpy.transpose(p))
    ncoords=numpy.array(newa)
    ncoords.shape=(-1,4)
    temp=numpy.take(ncoords,(0,1,2),-1)

    #qsel = numpy.arange(natoms-1)
    qsel = numpy.arange(natoms)
    totsnz1=numpy.zeros(len(temp)*3,numpy.int32)
    lts=0
    for ts in range(len(temp)):
        totsnz1[lts]=(qsel[ts]*3)
        totsnz1[lts+1]=(qsel[ts]*3)+1
        totsnz1[lts+2]=(qsel[ts]*3)+2
        lts=lts+3

    numpy.put(linker_coor[frame],totsnz1,temp)

    return	

def build_hybrid(linker,three_ball,minimum_base_angle):

    #### HARDWIRED TO BUILD FAB WITH A SET LINKER LENGTH!!!!
    #### HARDWIRED TO BUILD FAB WITH A SET LINKER LENGTH!!!!
    #### HARDWIRED TO BUILD FAB WITH A SET LINKER LENGTH!!!!
    #### OPEN	
    frame = 0

    hybrid = sasmol.SasHybrid(0)
    natoms = 2*linker.natoms() + three_ball.natoms()
    print 'hybrid.natoms() = ',natoms
    hybrid.setNatoms(natoms)

    coor = numpy.zeros((1,natoms,3),numpy.float32)

    get_pdb_values(hybrid,natoms)

    ### fc ball --> linker_z --> fab ball(z) 
    ### then
    ### --> rotated_linker --> fab ball(rotated)

    segname = hybrid.segname()
    name = hybrid.name() ; linker_name = linker.name()
    linker_resid = linker.resid()
    linker_resname = linker.resname()

    #segname[0]='FC' ; segname[16]='FAB1' ; segname[32]='FAB2'
    segname[0]='FAB1' ; segname[16]='FAB2' ; segname[32]='FC'
    count = 0
    for i in xrange(1,16):
        name[i] = linker_name[count]
        segname[i] = 'LNK1'	
        count += 1
    count = 0
    for i in xrange(17,32):
        name[i] = linker_name[count]
        segname[i] = 'LNK2'	
        count += 1
    hybrid.setSegname(segname)
    hybrid.setName(name)

    fc_coor = three_ball.coor()[0][0]
    fab1_coor = three_ball.coor()[0][1]
    fab2_coor = three_ball.coor()[0][2]

    mass = numpy.zeros(natoms,numpy.float32)

    linker_coor = linker.coor()[0]
    #coor[0][0] = fc_coor
    #coor[0][16] = fab1_coor
    #coor[0][32] = fab2_coor

    coor[0][0] = fab1_coor
    coor[0][16] = fab2_coor
    coor[0][32] = fc_coor

    for i in xrange(linker.natoms()):
        coor[0][i+1] = linker_coor[i]	
        mass[i+1] = 110.0
    linker.rotate(frame,'y',minimum_base_angle)
    for i in xrange(linker.natoms()):
        coor[0][i+17] = linker_coor[i]	
        mass[i+17] = 110.0

    hybrid.setCoor(coor)

    mass[0] = 49413.0
    mass[16] = 49413.0
    mass[32] = 49413.0
    hybrid.setMass(mass)
    mol_type = []
    for i in xrange(natoms):
        mol_type.append("protein")
    hybrid.setMoltype(mol_type)

    resid=numpy.zeros(natoms,numpy.int)
    resid[0]=1
    resname=[]
    resname.append("MAB")
    for i in xrange(1,16):
        resid[i] = linker_resid[i-1]+1
        resname.append(linker_resname[i-1])
        resid[i+16] = linker_resid[i-1]+7
    resid[16]=7 ; resname.append("MAB")
    resid[32]=13
    for i in xrange(1,16):
        resname.append(linker_resname[i-1])
    resname.append("MAB")

    hybrid.setResid(resid)
    hybrid.setResname(resname)

    hybrid.write_pdb('hybrid.pdb',frame,'w')

    return hybrid

def get_linker_assign_coor(three_ball,ball_diameter,minimum_base_angle,surface1_to_anchor_cord_length,linker_pdb):

    frame = 0 

    linker = sasmol.SasMol(0)
    linker.read_pdb(linker_pdb)

    # move linker com to origin

    linker.moveto(frame,[0.0,0.0,0.0])

    last_c = linker.coor()[0][-1]		

    # move linker last_c to origin

    linker.translate(frame,-last_c)  ### not debugged for all directions of linker ...

    first_n = linker.coor()[0][0] ; last_c = linker.coor()[0][-1]		
    linker_length = numpy.sqrt(numpy.sum((first_n-last_c)**2.0))

    # define end-point vector b

    b = first_n - last_c

    # define axis to rotate to (z-axis here)

    a = numpy.zeros(3,numpy.float)
    a[2] = 1.0

    # align coordinates onto axis (z-axis from "a" above)

    R = get_alignment_rotation_matrix(b,a)
    rotate_coordinates(b,R,linker)

    # move linker so that c-terminal is at origin

    linker.translate(frame,-linker.coor()[0][-1])

    # move linker to anchor point

    linker.moveto(frame,[0.0,0.0,(ball_diameter/2.0)+surface1_to_anchor_cord_length+(linker_length/2.0)])

    linker.write_pdb('new_single_linker.pdb',frame,'w')

    # build hybrid molecule

    hybrid = build_hybrid(linker,three_ball,minimum_base_angle)

    return hybrid

def angle_initialization(hybrid,rotation_variables):

    molecule_type = rotation_variables[0]
    temperature = rotation_variables[1]
    basis = rotation_variables[2]
    number_of_ranges = rotation_variables[3]
    reslow = rotation_variables[4]
    numcont = rotation_variables[5]
    dtheta = rotation_variables[6]


    ##### INITIALIZE VARIABLES FOR SEQUENCE THAT WILL BE ROTATED #####

    dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid = initialize_rotation.initialize_system(hybrid,molecule_type,basis,number_of_ranges,reslow,numcont)

    return dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid 


def main(simulation_type,number_of_steps,initial_pdb,output_pdb,argon_file_name,number_of_molecules,argon_box_length,length_scale_factor,ball_diameter,fc_to_fab_vector_length,minimum_base_angle, surface1_to_anchor_cord_length,linker_pdb,rotation_variables,number_of_groups,residues_in_groups,mab_diameter):

    start_type = 1
    coor_start_type = 2

    if(ball_diameter>0):
        cutoff = ball_diameter
    else:
        cutoff = 2.0

    if(start_type != 1):

    # allow for restart
        single_molecule = setup_mmc.get_base_molecule(initial_pdb,start_type)

    else:
        three_ball = build_ball_coordinates(ball_diameter,fc_to_fab_vector_length,minimum_base_angle)

        single_molecule = get_linker_assign_coor(three_ball,ball_diameter,minimum_base_angle,surface1_to_anchor_cord_length,linker_pdb)


    single_molecule._totalmass = 150000
    print 'hybrid mass (amu) = ',single_molecule.totalmass()

    dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid = angle_initialization(single_molecule,rotation_variables)

    groups,group_masks = group_rotation.make_groups(single_molecule,number_of_groups,residues_in_groups)

    if(simulation_type == 'single_mab'):

        simulate_single_mab(single_molecule,number_of_steps,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff)

        sys.exit()

    elif(simulation_type == 'box_of_mab'):

        all_molecules = sasmol.SasMol(1)

        if(start_type == 1):

            com_coor,boxl,max_diameter_molecule = get_com_coords_and_build_box(single_molecule,coor_start_type,number_of_molecules,argon_file_name,argon_box_length,length_scale_factor,mab_diameter)

            print 'len(com_coor) = ',len(com_coor)
            print com_coor

            #setup_mmc.duplicate_molecule_and_assign_coor(all_molecules,single_molecule,com_coor)
            duplicate_molecule_and_assign_coor(all_molecules,single_molecule,com_coor,initial_pdb,number_of_steps,groups,group_masks,number_of_groups,residues_in_groups,rotation_variables,dihedral_parameters,flexible_residues,residue_rotation_indices,residue_rotation_mask,step_parameters,first_last_resid,cutoff)

#			total_mass = all_molecules.calcmass() ; print ' > total mass = ',total_mass

#			all_molecules.write_pdb(output_pdb,0,'w')

if __name__ == '__main__':

    ##### REQUIRED INPUT AND SETUP #####

    simulation_type = 'single_mab'
    #simulation_type = 'box_of_mab'

    number_of_steps = 100

    #### VARIABLES FOR ROTATION

    molecule_type = 'protein'
    temperature = 300.0
    basis = 'C'
    number_of_ranges = 2 
    reslow = [3,9]
    numcont = [3,3]
    dtheta = [30.0,30.0]

    rotation_variables = [molecule_type,temperature,basis,number_of_ranges,reslow,numcont,dtheta]

    #### VARIABLES FOR GROUP DEFINITION

    number_of_groups = 2
    residues_in_groups = [[1,2,3,4,5,6], [7,8,9,10,11,12]]

    #### RESTART PDB

    initial_pdb = 'ball.pdb'

    #### LINKER PDB --> note that this code is set up so that both linkers are identical
    #### OPEN	

    linker_pdb = 'backbone_tverk.pdb'


    #### VARIABLES FOR BUILDING BOX OF MOLECULES: NOT NEEDED FOR SINGLE MAB

    output_pdb = '2048_balls.pdb'
    argon_file_name = 't_85_p_419_nvt4.coor'
    number_of_molecules = 2048
    argon_box_length = 45.0131
    length_scale_factor = 2.0/0.21
    diameter_argon = 1.88*2.0
    length_scale_factor = 70.0/diameter_argon

    mab_diameter = 90.0  # this is the default minimum diameter

    ### INFORMATION FOR BUILDING MAB SPECIFIC HYBRID

    ###
    ### from the center of one ball to the center of the other ball 
    ### from the center of Fc to center of Fab
    ###
    ### ball_radius + surface1_to_anchor_cord_length + maximum_anchor_to_surface2_cord_length + ball_radius
    ### 
    ### where the angle between the two vectors is given by minimum_base_angle
    ###

    diameter_carbon = 3.4
    ball_diameter = 16.0*diameter_carbon

    ball_radius = ball_diameter/2.0
    surface1_to_anchor_cord_length = 18.0
    maximum_anchor_to_surface2_cord_length = 14.46

    fc_to_fab_vector_length = ball_radius + surface1_to_anchor_cord_length + maximum_anchor_to_surface2_cord_length + ball_radius

    minimum_base_angle = 36.5 * numpy.pi/180.0

    print 'fc_to_fab_vector_length = ',fc_to_fab_vector_length
    print 'minimum_base_angle ',minimum_base_angle


    ### MAIN PROGRAM CALL

    main(simulation_type,number_of_steps,initial_pdb,output_pdb,argon_file_name,number_of_molecules,argon_box_length,length_scale_factor,ball_diameter,fc_to_fab_vector_length,minimum_base_angle, surface1_to_anchor_cord_length,linker_pdb,rotation_variables,number_of_groups,residues_in_groups,mab_diameter)


    ### JUNK

    '''
	volume_of_ball = 4.0*numpy.pi*((ball_diameter/2.0)**3)/3.0
	protein_density_g_per_cm3 = 1.47*0.7
	cm3_to_angstrom3 = 1E-24
	g_to_amu = 1.0/(1000*1.660538E-27)
	protein_density = protein_density_g_per_cm3 * cm3_to_angstrom3 * g_to_amu
	linker_mass = 16 * 110.0
	single_molecule._totalmass = protein_density*3*volume_of_ball + linker_mass
	'''