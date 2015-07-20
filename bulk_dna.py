#!/usr/bin/python
# $Id$
# time using FORTRAN double loop, N=1000, iters=1000 (so 1*10^6 steps): 958.887075186 seconds
# time using python double loop, N=1000, iters=1000 (so 1*10^6 steps):

import sassie.sasmol.sasmol as sasmol
import sassie.simulate.monte_carlo.double_stranded_nucleic as double_stranded_nucleic
import sassie.util.module_utilities as module_utilities
import multiprocessing
from sassie.simulate.monte_carlo.double_stranded_nucleic import *
import sassie.interface.input_filter as input_filter
import sassie.simulate.monte_carlo.monte_carlo as monte_carlo

import numpy as np
import time
# improt string, os, locale, sys, random

DEBUG = True
app = 'long_ds_dna_test'


def rotate_coarse_grained_beads(other_self, group_number):
    '''
    method to choose a cg bead to rotate
    '''

    log = other_self.log
    log.debug('in rotate_coarse_grained_beads')

    vars = other_self.vars
    nvars = other_self.vars.nvars
    pgui = other_self.run_utils.print_gui
    group_mol = other_self.group_molecules[group_number]
    run_path = other_self.runpath

    n_flex_regions = vars.number_of_flexible_regions
    theta_max = vars.delta_theta_array[group_number]

    cg_flex_mol = nvars.coarse_grained_flexible_mol[group_number]
    post_mol = nvars.post_mol[group_number]
    xyz_vectors = nvars.xyz_vectors[group_number]
    all_beads = nvars.all_beads[group_number]
    flexible_mol = nvars.flexible_mol[group_number]
    all_overlap_mask = nvars.all_overlap_mask
    post_overlap_mask = nvars.post_overlap_mask[group_number]
    bead_masks = nvars.nucleic_acid_bead_masks[group_number]
    post_mask = nvars.post_mask[group_number]
    flexible_mask = nvars.flexible_mask[group_number]

    try:
        soft_rotation = nvars.soft_rotation
    except:
        soft_rotation = nvars.soft_rotation = 1
    if soft_rotation > 1:
        log.info('softening all rotations over %d coarse-grained beads' %
                 soft_rotation)
    try:
        i_loop = nvars.i_loop
    except:
        i_loop = nvars.i_loop = 1

    try:
        persistence_length = nvars.persistence_length
    except:
        persistence_length = nvars.persistence_length = 530.0  # 530 A

    try:
        ds_na_type = nvars.double_stranded_nucleic_acid_type
    except:
        ds_na_type = nvars.double_stranded_nucleic_acid_type = 'b'

    try:
        theta_z_max = nvars.delta_theta_z_array[group_number]
    except:
        nvars.delta_theta_z_array = [1.0] * vars.number_of_flexible_regions
        theta_z_max = nvars.delta_theta_z_array[group_number]

    try:
        trials = nvars.trials
    except:
        trials = nvars.trials = 1

    try:
        dcd_write_frequency = nvars.dcd_write_frequency
    except:
        dcd_write_frequency = nvars.dcd_write_frequency = 1

    try:
        n_go_back = nvars.n_go_back
    except:
        n_go_back = nvars.n_go_back = 50

    try:
        keep_unique = nvars.keep_unique
    except:
        keep_unique = nvars.keep_unique = True

    try:
        keep_cg_files = nvars.keep_cg_files
    except:
        keep_cg_files = nvars.keep_cg_files = True

    try:
        min_frequency = nvars.min_frequency
    except:
        min_frequency = nvars.min_frequency = 100

    n_cg_beads = cg_flex_mol.natoms()

    ### not needed or long moves ###
    # s group_all_atom_out_file = '%s/na_group%d_%03d.dcd' % (run_path,
    # s                                                      group_number, i_loop)
    # s log.info('opening new dcd file to store trajectory: %s' % group_all_atom_out_file)
    # s group_all_atom_dcd_out =
    # group_mol.open_dcd_write(group_all_atom_out_file)

    # s # create the coarse-grained nucleic acid and post dcd and pdb files
    # s cg_flex_out_file = '%s/cg_na_group%d_%03d.dcd' % (run_path, group_number,
    # s                                                  i_loop)
    # s log.info('opening new pdb file to store original go_back structure: '
    # s          '%s.pdb' % cg_flex_out_file[:-4])
    # s cg_flex_mol.write_pdb(cg_flex_out_file[:-4] + '.pdb', 0, 'w')
    # s log.info('opening new dcd file to store trajectory: %s' % cg_flex_out_file)
    # s cg_flex_dcd_out = cg_flex_mol.open_dcd_write(cg_flex_out_file)

    # s post_out_file = '%s/post_group%d_%03d.dcd' % (run_path, group_number,
    # s                                              i_loop)
    # s log.info('opening new pdb file to store original go_back structure: '
    # s          '%s.pdb' % post_out_file[:-4])
    # s post_mol.write_pdb(post_out_file[:-4] + '.pdb', 0, 'w')
    # s log.info('opening new dcd file to store trajectory: %s' % post_out_file)
    # s post_dcd_out = post_mol.open_dcd_write(post_out_file)
    # s
    # s # will write these out to dcd files to store the coordinates along the
    # way
    x_vec_mol = sasmol.SasMol(0)
    y_vec_mol = sasmol.SasMol(0)
    z_vec_mol = sasmol.SasMol(0)
    # s error, mask = cg_flex_mol.get_subset_mask('(all)')
    # s handle_error(other_self, error)
    # s error = cg_flex_mol.copy_molecule_using_mask(x_vec_mol, mask, 0)
    # s handle_error(other_self, error)
    # s error = cg_flex_mol.copy_molecule_using_mask(y_vec_mol, mask, 0)
    # s handle_error(other_self, error)
    # s error = cg_flex_mol.copy_molecule_using_mask(z_vec_mol, mask, 0)
    # s handle_error(other_self, error)
    # the np.array recast these so they
    x_vec_mol.setCoor(np.array([xyz_vectors[0]]))
    y_vec_mol.setCoor(np.array([xyz_vectors[1]]))  # do not update with vecXYZ
    z_vec_mol.setCoor(np.array([xyz_vectors[2]]))
    # s x_vec_dcd_file = '%s/x_vec_%03d.dcd' % (run_path, i_loop)
    # s y_vec_dcd_file = '%s/y_vec_%03d.dcd' % (run_path, i_loop)
    # s z_vec_dcd_file = '%s/z_vec_%03d.dcd' % (run_path, i_loop)
    # s log.info('opening new dcd file to store coordinate vectors: %s' % x_vec_dcd_file)
    # s x_vec_dcd_out = x_vec_mol.open_dcd_write(x_vec_dcd_file)
    # s log.info('opening new dcd file to store coordinate vectors: %s' % y_vec_dcd_file)
    # s y_vec_dcd_out = y_vec_mol.open_dcd_write(y_vec_dcd_file)
    # s log.info('opening new dcd file to store coordinate vectors: %s' % z_vec_dcd_file)
    # s z_vec_dcd_out = z_vec_mol.open_dcd_write(z_vec_dcd_file)
    # s x_vec_mol.write_dcd_step(x_vec_dcd_out, 0, 1)
    # s y_vec_mol.write_dcd_step(y_vec_dcd_out, 0, 1)
    # s z_vec_mol.write_dcd_step(z_vec_dcd_out, 0, 1)

    # initialize variables for each run
    steps_from_0 = np.zeros(trials, dtype='int64')

    xyz = np.copy(xyz_vectors)
    flex_coor = np.copy(cg_flex_mol.coor()[0])  # unique memory for each

    post_coor = np.copy(post_mol.coor()[0])  # unique memory for each

    (bead_unit_vectors, mean_bead_distance) = get_chain_parameters(other_self,
                                                                   flex_coor)

    bead_radius = np.floor(mean_bead_distance.min())
    post_radius = 0.4                                 # 1 A

    # effective width of the ds nucleic acid chain
    # 4.6nm for B-Form (taken from: Odijk, T. Phys. Rev. E 2008, 77, 060901(R))
    # yet to use a, and z type nucleic acid
    na_energetics_width = {'a': 0, 'b': 46., 'z': 0}

    try:
        nvars.chain_width[group_number]
    except:
        nvars.chain_width = [None] * n_flex_regions
    chain_width = nvars.chain_width[group_number]
    if chain_width is None:
        chain_width = na_energetics_width[ds_na_type.lower()]
        if chain_width > mean_bead_distance.any():
            chain_width = bead_radius
            log.debug('chain width set to %d Angstroms (so it is less than the '
                      'distance between beads)' % chain_width)
        nvars.chain_width[group_number] = chain_width

    # only need to consider collision between stationary NA and post molecules
    bead_post_test = bead_radius + post_radius

    # calculate the energy of the starting positions
    wca_energy_matrix_old = np.zeros((n_cg_beads, n_cg_beads))
    bend_energy_old = bend_energy(
        other_self, persistence_length, bead_unit_vectors, mean_bead_distance)

    (wca_energy_old, wca_energy_matrix_old) = get_wca_energy(other_self,
                                                             chain_width, flex_coor, wca_energy_matrix_old, 0)

    total_energy_old = bend_energy_old + wca_energy_old

    n_accept = 0  # total times configuration was accepted
    n_reject = 0  # total times configuration was rejected
    n_written = 0  # total times dcd write has been called
    fail_tally = 0  # number of times failed for particular iteration
    n_from_reload = 0  # number of stps since last reload
    n_reload = [0]  # list containing the i_goback values

    # s while n_accept < trials: # for calculating the Rg and Re we want every
    # fails counted as separate structures
    for trial in xrange(trials):
        # Choose a bead to rotate
        # rotating bead 0 doesn't mean anything
        trial_bead = np.random.randint(1, n_cg_beads)

        # Determine rotation to perform
        theta_max_xyz = np.array(
            [theta_max, theta_max, theta_z_max])  # degrees
        theta_xyz = np.zeros(3)
        # theta_xyz[np.random.randint(3)] = 1.0 # randomly select x, y, or z
        # randomly select x, or y   # do not want to use z rotations
        theta_xyz[np.random.randint(2)] = 1.0
        # continuously distribution
        theta_xyz *= theta_max_xyz * (2 * np.random.random() - 1)
        theta_xyz /= soft_rotation
        log.debug('[x, y, z] rotation angle: [%f, %f, %f]' %
                  (theta_xyz[0], theta_xyz[1], theta_xyz[2]))

        # generate a newly rotated model
        #!# may need to modify if there is no post coordinates #TEST#
        (flex_coor[trial_bead:], xyz[:, trial_bead:], post_coor) = bead_rotate(
            other_self, flex_coor[trial_bead - 1:], xyz[:, trial_bead - 1:],
            theta_xyz, post_coor, soft_rotation)

        # calculate the change in energy (dU) and the boltzman factor (p)
        (bead_unit_vectors, bead_distances) = get_chain_parameters(
            other_self, flex_coor)
        bend_energy_new = bend_energy(other_self, persistence_length,
                                      bead_unit_vectors, bead_distances)

        # ~~~~ DNA interaction energy  ~~~~~~#
        (wca_energy_new, wca_energy_matrix_new) = get_wca_energy(other_self,
                                                                 chain_width, flex_coor, wca_energy_matrix_old, trial_bead)

        total_energy_new = bend_energy_new + wca_energy_new
        energy_change = total_energy_new - total_energy_old

        with warnings.catch_warnings():
            warnings.filterwarnings('error')  # need this for np warnings
            try:
                probability = np.exp(-energy_change)
                log.debug('probability of new structure: %f' % probability)
            except Warning:
                if energy_change > 99:
                    probability = 0
                    log.debug('large energy increase, probability set to 0')
                elif energy_change < 0:
                    probability = 1
                    log.debug('energy decreased, probability set to 1')
            except:
                log.error("ERROR: Unexpected error calculating Boltzmann "
                          "probability:" + sys.exc_info()[0])

        test = np.random.random()
        collision = 0

        if test > probability:  # pass: test <= probability
            # if False:
            ds_nucleic_acid_pass = False
            log.debug(
                'step rejected: ds nucleic acid bending energy parameters')
        else:
            ds_nucleic_acid_pass = True

            error, post_heavy_coor = post_mol.get_coor_using_mask(
                0, post_overlap_mask)
            handle_error(other_self, error)
            post_heavy_coor = post_heavy_coor[0]
            # check the group for overlap
            # only need to consider collision between stationary NA and post
            # stationary NA:      flex_coor[trial_bead:]
            # post (heavy atoms): post_heavy_coor
            if 1 == two_body_overlap(other_self, flex_coor[trial_bead:],
                                     post_heavy_coor, bead_post_test):
                log.debug(
                    'step rejected: collision between post atoms and the fixed NA beads')
                collision = 1

        if ds_nucleic_acid_pass and collision == 0:
            n_from_reload += 1
            steps_from_0[n_accept] = n_from_reload + n_reload[-1]
            # increment accept counter
            n_accept += 1
            # update flex coordinates
            cg_flex_mol.setCoor(np.array([flex_coor]))
            # update post coordinates
            post_mol.setCoor(np.array([post_coor]))
            # increment the number of moves for the group
            nvars.n_moves[group_number] += 1
            xyz_vectors = np.copy(xyz)                 # update NA orientations
            x_vec_mol.setCoor(np.array([xyz_vectors[0]]))
            y_vec_mol.setCoor(np.array([xyz_vectors[1]]))
            z_vec_mol.setCoor(np.array([xyz_vectors[2]]))

            wca_energy_matrix_old = np.copy(
                wca_energy_matrix_new)  # update dsNA WCA energy
            # update total energy
            total_energy_old = total_energy_new

            log.debug('trial_bead(%3d) = %2d\t failed attempts = %2d' % (
                n_accept, trial_bead, fail_tally))
            log.info('.')
            fail_tally = 0                     # reset fail_tally

            # recover an all atom representation and save coordinates to a dcd
            if 0 == n_accept % dcd_write_frequency:
                # ~~recover aa-DNA~~
                recover_all_atom_dsNA(other_self, cg_flex_mol, flexible_mol,
                                      xyz_vectors, all_beads, bead_masks)

                # ~~recover aa-Protein~~
                # recover_aaPro_model(aa_pgroup_masks, rigid_group_masks, rigid_mol,
                # rigid_group_mols, aa_pro)

                # ~~Combine aa Complete Structure~~
                group_mol.set_coor_using_mask(post_mol, 0, post_mask)
                group_mol.set_coor_using_mask(flexible_mol, 0, flexible_mask)
                # ~~Write DCD step~~
                n_written += 1
                group_mol.write_dcd_step(group_all_atom_dcd_out, 0, n_written)

                # write out the accepted configuration for go-back use
                if n_go_back > 0:
                    # to turn off, n_go_back should be set to negative value

                    cg_flex_mol.write_dcd_step(cg_flex_dcd_out, 0, n_written)
                    post_mol.write_dcd_step(post_dcd_out, 0, n_written)
                    # these are incremented by one because the 0th contains the
                    # original coordinates
                    x_vec_mol.write_dcd_step(x_vec_dcd_out, 0, n_written + 1)
                    y_vec_mol.write_dcd_step(y_vec_dcd_out, 0, n_written + 1)
                    z_vec_mol.write_dcd_step(z_vec_dcd_out, 0, n_written + 1)
        else:
            if fail_tally == n_go_back:
                # load a previous state and return the state index
                i_go_back = go_back(other_self, n_accept, cg_flex_out_file,
                                    cg_flex_mol, post_out_file, post_mol,
                                    x_vec_dcd_file, x_vec_mol, y_vec_mol,
                                    y_vec_dcd_file, z_vec_mol, z_vec_dcd_file,
                                    xyz_vectors)
                # get the flex_coor from the loaded file
                flex_coor = np.copy(cg_flex_mol.coor()[0])

                # get the reference energy for the loaded configuration
                (bead_unit_vectors, bead_distances) = get_chain_parameters(
                    other_self, flex_coor)
                (wca_energy_old, wca_energy_matrix_old) = get_wca_energy(
                    other_self, chain_width, flex_coor, wca_energy_matrix_old, 0)
                bend_energy_old = bend_energy(
                    other_self, persistence_length, bead_unit_vectors, bead_distances)
                total_energy_old = bend_energy_old + wca_energy_old

                n_from_reload = 0
                n_reload.append(steps_from_0[i_goback - 1])
                fail_tally = 0  # reset the fail counter
            else:
                fail_tally += 1                 # increment bead reject counter
                # increment total reject counter
                n_reject += 1
                # reset the flexible coordinates
                flex_coor = np.copy(cg_flex_mol.coor()[0])

            # reset the post coordinates
            post_coor = np.copy(post_mol.coor()[0])
            # reset the bead orientations
            xyz = np.copy(xyz_vectors)

            # save previous coordinates again
            if not keep_unique:
                # ~~Write DCD step~~
                n_written += 1
                group_mol.write_dcd_step(group_all_atom_dcd_out, 0, n_written)

    # s group_mol.close_dcd_write(group_all_atom_dcd_out)

    if 0 == nvars.n_moves[group_number] % min_frequency:
        energy_minimize_nucleic(other_self, group_number)
    else:
        # Need to decide on when to do energy minimization
        pass

    try:
        log.info(
            'removing dcd file used to store coordinate vectors: %s' % x_vec_dcd_file)
        os.remove(x_vec_dcd_file)
        log.info(
            'removing dcd file used to store coordinate vectors: %s' % y_vec_dcd_file)
        os.remove(y_vec_dcd_file)
        log.info(
            'removing dcd file used to store coordinate vectors: %s' % z_vec_dcd_file)
        os.remove(z_vec_dcd_file)

        if keep_cg_files:
            cg_flex_mol.close_dcd_write(cg_flex_dcd_out)
            post_mol.close_dcd_write(post_dcd_out)
        else:
            log.info(
                'removing coarse-grained flexible nucleic acid pdb file stored for go_back: %s.pdb' % cg_flex_out_file[:-4])
            os.remove(cg_flex_out_file[:-4] + '.pdb')
            log.info(
                'removing coarse-grained flexible nucleic acid dcd file stored for go_back: %s' % cg_flex_out_file)
            os.remove(cg_flex_out_file)
            log.info(
                'removing post pdb file stored for go_back: %s.pdb' % post_out_file[:-4])
            os.remove(post_out_file[:-4] + '.pdb')
            log.info('removing post dcd file stored for go_back: %s' %
                     post_out_file)
            os.remove(post_out_file)
    except:
        pass

    if n_go_back > 0 and DEBUG:
        go_back_out_file = '%s/n_from_0_%03d.txt' % (run_path, i_loop)
        log.debug('saving go_back output to %s file' % go_back_out_file)
        np.savetxt(go_back_out_file, steps_from_0, fmt='%d')

    log.debug('accepted %d moves' % n_accept)
    log.debug('rejected %d moves' % n_reject)

    log.debug('accepted %0.2f percent of trials' %
              (100.0 * n_accept / (n_reject + n_accept)))
    nvars.n_accept = n_accept
    nvars.n_reject = n_reject
    return


def makeLongDNA(n_lp):
    print 'making DNA that is %d*lp long' % n_lp

    # 15 bp/bead or 51 A/bead (3.4 A/bp)

    lp = 530  # persistence length in A
    l = 2**(1. / 6.) * 46   # separation distance between beads = 51.6A

    longDNA = sasmol.SasMol(0)
    L = n_lp * lp
    N = int(L / l)
    natoms = N + 1
    print 'natoms = ', natoms
    longDNA._L = L
    longDNA._natoms = natoms

    # initialize the long DNA coordinates
    longCoor = np.zeros((1, natoms, 3), np.float)
    # set the z-values to the index of the array
    longCoor[0][:, 2] = range(natoms)
    # scale the z-values to the right seperation
    longCoor *= l
    # print longCoor[-5:]

    longDNA.setCoor(longCoor)
    longDNA.setElement(['C'] * natoms)

    xyz_vectors = np.zeros((3, natoms, 3))
    xyz_vectors[0] = [1, 0, 0]
    xyz_vectors[1] = [0, 1, 0]
    xyz_vectors[2] = [0, 0, 1]

    # n = L/l                         # number of times need to repeat the grain
    # print '(l, L, n)', (l, L, n)

    return (longDNA, xyz_vectors)


def unpack_variables(other_self, variables):
    '''
    method to extract variables into system wide class instance
    '''

    other_self.log.debug('in unpack_variables')

    vars.runname = variables['runname'][0]
    vars.dcdfile = variables['dcdfile'][0]
    vars.pdbfile = variables['pdbfile'][0]
    vars.number_of_flexible_regions = variables[
        'number_of_flexible_regions'][0]
    vars.basis_string_array = variables['basis_string_array'][0]
    vars.delta_theta_array = variables['delta_theta_array'][0]
    vars.rotation_type_array = variables['rotation_type_array'][0]
    vars.rotation_direction_array = variables['rotation_direction_array'][0]
    vars.basis_overlap_array = variables['basis_overlap_array'][0]
    vars.post_basis_string_array = variables['post_basis_string_array'][0]
    vars.post_basis_overlap_array = variables['post_basis_overlap_array'][0]
    vars.temperature = variables['temperature'][0]
    vars.trial_steps = variables['trial_steps'][0]
    vars.goback = variables['goback'][0]
    vars.cutoff = variables['cutoff'][0]

    vars.nonbondflag = variables['nonbondflag'][0]
    vars.seed = variables['seed'][0]

    self.vars = vars

    return


def mag(vec):

    (r, ) = vec.shape
    sumSquares = 0.0
    for i in xrange(r):
        sumSquares += vec[i]**2

    return np.sqrt(sumSquares)


class vars():

    def __init__(self, parent=None):
        pass


class nvars():

    def __init__(self, parent=None):
        pass

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
    nvars.trials = 10   # intend to use 1E5
    delta_theta_array = '15'
    Llp = 15    # L/lp
    # ----- Modify these ---

    # ----- DO NOT Modify these ---
    nvars.persistence_length = 530      # persistence length  (lp = 530A)
    nvars.width = [46]        # width of chain in A (w = 46A)l
    # ----- DO NOT Modify these ---
    lp = nvars.persistence_length

    svariables = {}
    svariables['runname'] = ('%dlp_dnaMoves' % Llp, 'string')
    svariables['runname'] = ('', 'string')
    svariables['number_of_flexible_regions'] = ('1', 'int')
    svariables['temperature'] = ('300.0', 'float')
    svariables['seed'] = ('0, 123', 'int_array')
    svariables['delta_theta_array'] = (delta_theta_array, 'float_array')

    svariables['dcdfile'] = ('dummy', 'string')
    svariables['pdbfile'] = ('dummy', 'string')
    svariables['basis_string_array'] = ('dummy', 'string')
    svariables['rotation_type_array'] = ('dummy', 'string')
    svariables['rotation_direction_array'] = ('dummy', 'string')
    svariables['basis_overlap_array'] = ('dummy', 'string')
    svariables['post_basis_string_array'] = ('dummy', 'string')
    svariables['post_basis_overlap_array'] = ('dummy', 'string')
    svariables['trial_steps'] = ('-1', 'int')
    svariables['goback'] = ('-1', 'int')
    svariables['cutoff'] = ('4.0', 'float')
    svariables['nonbondflag'] = ('0', 'int')

    error, variables = input_filter.type_check_and_convert(svariables)

    self = monte_carlo.simulation()
    self.txtOutput = multiprocessing.JoinableQueue()

    self.run_utils = module_utilities.run_utils(app, self.txtOutput)
    self.run_utils.setup_logging(self)
    unpack_variables(self, variables)
    self.run_utils.general_setup(self)
    self.group_molecules = [sasmol.SasMol(0)]
    log = self.log

    dummy = sasmol.SasMol(0)
    dummy.read_pdb('dummy.pdb')
    never = nvars.trials + 1
    nvars.n_flex_regions = 1
    nvars.soft_rotation = 1
    nvars.i_loop = 0
    nvars.dcd_write_frequency = never
    nvars.n_go_back = -1
    nvars.keep_cg_files = False
    nvars.min_frequency = never
    nvars.post_mol = [dummy]
    nvars.flexible_mol = [dummy]
    nvars.post_mol = [dummy]
    nvars.all_beads = [dummy]
    nvars.all_overlap_mask = [np.array([0])]
    nvars.post_overlap_mask = [np.array([0])]
    nvars.nucleic_acid_bead_masks = [np.array([0])]
    nvars.post_mask = [np.array([0])]
    nvars.flexible_mask = [np.array([0])]
    nvars.n_moves = [0]
    self.group_molecules = [dummy]

    (cg_dna, xyz_vectors) = makeLongDNA(Llp)
    nvars.coarse_grained_flexible_mol = [cg_dna]
    nvars.xyz_vectors = [xyz_vectors]
    natoms = cg_dna.natoms()
    cg_dna.setAtom(['ATOM'] * natoms)
    cg_dna.setIndex(np.array(range(1, natoms + 1)))
    cg_dna.setResid(np.array(range(1, natoms + 1)))
    cg_dna.setLoc([' '] * natoms)
    cg_dna.setResname(['ADE'] * natoms)
    cg_dna.setChain(['A'] * natoms)
    cg_dna.setRescode([' '] * natoms)
    cg_dna.setOccupancy(['0.00'] * natoms)
    cg_dna.setBeta(['0.00'] * natoms)
    cg_dna.setSegname(['DNA'] * natoms)
    cg_dna.setElement(['C'] * natoms)
    cg_dna.setName(['C'] * natoms)
    cg_dna.setCharge(['  '] * natoms)

    rg0 = cg_dna.calcrg(0) / lp
    h0 = mag(cg_dna.coor()[0, -1] - cg_dna.coor()[0, 0]) / lp

    log.info('iter: %d of %d; (a, r) = (  ,  ); rg/lp = %0.2f; h/lp = %0.2f' %
             (0, iters, rg0, h0))

    timestr = time.strftime("%y%m%d_%H%M%S")
    fName = self.runpath + '/' + timestr + '_%dlp.o' % Llp

    cg_dna.write_pdb(fName[:-1] + 'pdb', 0, 'w')
    cg_dna.open_dcd_write(fName[:-1] + 'dcd')
    outData = open(fName, 'a')
    log.info("# L=%d\t iters=%d\t nsteps=%d\t theta_max=%f \n# moves\t rg/lp\t\t re/lp\t\t a\t r\n" %
             (Llp, iters, nvars.trials, vars.delta_theta_array[0]))
    outData.close()

    tic = time.time()

    self.vars.nvars = nvars

    for i in xrange(iters):
        rotate_coarse_grained_beads(self, 0)

        rg_lp = cg_dna.calcrg(0) / lp
        h_lp = mag(cg_dna.coor()[0, -1] - cg_dna.coor()[0, 0]) / lp
        log.info('iter: %d of %d; (a, r) = (%d, %d); rg/lp = %0.2f; h/lp = %0.2f' %
                 (0, iters, nvars.n_accept, nvars.n_reject, rg_lp, h_lp))

        outData = open(fName, 'a')
        outData.write("%d\t%f\t%f\t%d\t%r\n" % (
            (i + 1) * nvars.trials, rg_lp, h_lp, nvars.n_accept, nvars.n_reject))
        outData.close()

        cg_dna.write_dcd_step(fName[:-1] + 'dcd', 0, i)

    cg_dna.close_dcd_write(fName[:-1] + 'dcd')
    toc = time.time() - tic
    log.debug('run time =', toc, 'seconds')

    # recover an all atom representation
    recover_aaDNA_model(cg_dna, xyz_vectors, allBeads)


'''
to do: 
reproduce one of the figures from D. Tree's paper (to verify the model is programmed the way they did)
if the model is fine and good for the length scale of interest then go to the mapping
plotting the energies 9


sassie: mixed monte carlo -> molecular
'''
# all_atom_pdb = 'trimer_stacked.pdb'

# s figure()
# s x = range(iters)
# s plot(x,rg)
# s xlabel('iteration #')
# s ylabel('Rg')

#test = np.eye(4)
# s est = np.zeros((10,4))
# s est[:,1] = range(10)
# s est[:,3] = 1
# s est[:,0] = range(10)
# s est[:,0] *= -1
# s or i in xrange(4):
# s        (T, Ti) = move2origin(test[i:])
# s        test[i:] = np.dot(test[i:],T)
# s        print 'T:\n', T
# s        print 'moved2origin:\n', test[i:]
# s
# s        (A, Ai) = align2z(test[i:])
# s        print 'before align: \n', test[i:]
# s        test[i:] = np.dot(test[i:],A)
# s        print 'x A ->\n', A
# s        print 'aligned: \n', test[i:]
# s        test[i:] = np.dot(test[i:],Ai)
# s        print 'un-aligned:  \n', test[i:]
# s
# s        test[i:] = np.dot(test[i:],Ti)
# s        print 'back2original:\n', test[i:]
# s
# s rint test
