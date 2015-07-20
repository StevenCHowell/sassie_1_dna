#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: drive the DNA moves toward a smaller (or larger) Rg
# Created: 9 October 2014
#
# $Id$
#
# 0000000011111111112222222222333333333344444444445555555555666666666677777777778
# 2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sassie.sasmol.sasmol as sasmol
import numpy as np
import dna.cgDNA_move as dna_move
import warnings
import time
import os


def driven_dna_mc(ARGS, cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, lp, trialbeads,
                  beadgroups, group_masks, all_beads, dna_bead_masks, aa_pgroup_masks,
                  cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask,
                  dna_type='b'):
    '''
    this function perform nsteps Monte-Carlo moves on the cg_dna
    '''

    timestr = time.strftime("%y%m%d_%H%M%S_")  # prefix for output files
    all_dcd_name = timestr + ARGS.pdb[:-4] + '.dcd'
    aa_all_dcd_out = aa_all.open_dcd_write(all_dcd_name)

    if False:
        aa_all.send_coordinates_to_vmd(2222, 0)

    # create the coarse-grained DNA and protein dcd and pdb files
    cg_dna_dcd_name = timestr + 'cg_dna.dcd'
    cg_pro_dcd_name = timestr + 'cg_pro.dcd'
    cg_dna_dcd_out = cg_dna.open_dcd_write(cg_dna_dcd_name)
    cg_pro_dcd_out = cg_pro.open_dcd_write(cg_pro_dcd_name)
    cg_dna.write_dcd_step(cg_dna_dcd_out, 0, 1)
    cg_pro.write_dcd_step(cg_pro_dcd_out, 0, 1)
    cg_dna.write_pdb(timestr + 'cg_dna.pdb', 0, 'w')
    cg_pro.write_pdb(timestr + 'cg_pro.pdb', 0, 'w')

    # create a dummy sasmol object for the 3 orientation vectors for each bead
    # will write these out to dcd files to store the coordinates along the way
    vecX_mol = sasmol.SasMol(0)
    vecY_mol = sasmol.SasMol(0)
    vecZ_mol = sasmol.SasMol(0)
    error, mask = cg_dna.get_subset_mask('(all)')
    error = cg_dna.copy_molecule_using_mask(vecX_mol, mask, 0)
    error = cg_dna.copy_molecule_using_mask(vecY_mol, mask, 0)
    error = cg_dna.copy_molecule_using_mask(vecZ_mol, mask, 0)
    # the np.array recast these so they
    vecX_mol.setCoor(np.array([vecXYZ[0]]))
    vecY_mol.setCoor(np.array([vecXYZ[1]]))  # do not update with vecXYZ
    vecZ_mol.setCoor(np.array([vecXYZ[2]]))
    vecX_dcd_name = timestr + 'vecX.dcd'
    vecY_dcd_name = timestr + 'vecY.dcd'
    vecZ_dcd_name = timestr + 'vecZ.dcd'
    vecX_dcd_out = vecX_mol.open_dcd_write(vecX_dcd_name)
    vecY_dcd_out = vecY_mol.open_dcd_write(vecY_dcd_name)
    vecZ_dcd_out = vecZ_mol.open_dcd_write(vecZ_dcd_name)
    vecX_mol.write_dcd_step(vecX_dcd_out, 0, 1)
    vecY_mol.write_dcd_step(vecY_dcd_out, 0, 1)
    vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, 1)

    # initialize variables for each run
    steps_from_0 = np.zeros(ARGS.nsteps, dtype='int64')
    xyz = np.copy(vecXYZ)
    d_coor = np.copy(cg_dna.coor()[0])  # unique memory for each
    p_coor = np.copy(cg_pro.coor()[0])  # unique memory for each

    # vectors between beads u, and average distance l
    (u, l) = dna_move.checkU(d_coor)
    # s print "(u, l) =", (u, l) # debug info
    lpl = lp / l  # setup the presistence length paramater

    # yet to use a, and z type dna
    dna_energy_width = {'a': 0, 'b': 46., 'z': 0}
    w = dna_energy_width[dna_type.lower()]
    if w > l:
        w = np.floor(l)
        # print '~~~ %.2f > %.2f ~~~~~~~~~~~~~~~~~~~~~~~~' % (w, l)
        print ('>>> setting chain width (w) to %d (chain width < distance' % w,
               ' btwn beads)')

    dna_diam = {'a': 25.5, 'b': 23.7, 'z': 18.4}
    dna_bead_radius = 4.5

    pro_bead_radius = 1.0  # 2A min seperation of CA atoms in database

    pro_pro_test = pro_bead_radius + pro_bead_radius
    dna_pro_test = dna_bead_radius + pro_bead_radius

    # calculate the energy of the starting positions
    wca0 = np.zeros((cg_dna.natoms(), cg_dna.natoms()))
    Ub0 = dna_move.energyBend(lpl, u, l)

    (Uwca0, wca0) = dna_move.f_energy_wca(w, d_coor, wca0, 0)

    U_T0 = Ub0 + Uwca0
    # print '(Ub0, Uwca0, Ub0/U_T0, Uwca0/U_T0) = ', (Ub0, Uwca0, Ub0/U_T0,
    #                                                 Uwca0/U_T0)
    n_accept = 0  # total times configuration was accepted
    n_reject = 0  # total times configuration was rejected
    n_written = 0  # total times dcd write has been called
    fail_tally = 0  # number of times failed for particular iteration
    n_from_reload = 0  # number of stps since last reload
    n_reload = [0]  # listt containing the i_goback values

    # this should not actually be >=, come back to this
    assert np.size(ARGS.theta_max) - 1 >= np.max(beadgroups), (
        'each group needs its own theta_max: %d < %d'
        % (np.size(ARGS.theta_max) - 1, np.max(beadgroups)))

    rg_old = cg_dna.calcrg(0)

    # Main MC loop #
    while n_accept < ARGS.nsteps:

        # Choose a bead to rotate
        trial_bead = trialbeads[int((trialbeads.size) * np.random.random())]

        # Determine rotation to perform
        theta_max = ARGS.theta_max[beadgroups[trial_bead]]
        # option to scale thetaZ separatly
        thetaZ_max = 0 * np.float(theta_max)
        # thetaZ_max = np.float(theta_max) # option to scale thetaZ separatly
        thetaZ = 2 * thetaZ_max * np.random.random() - thetaZ_max
        thetaX = 2 * theta_max * np.random.random() - theta_max
        thetaY = 2 * theta_max * np.random.random() - theta_max
        thetaXYZ = [
            thetaX / ARGS.n_soft, thetaY / ARGS.n_soft, thetaZ / ARGS.n_soft]
        # print theta_max, thetaXYZ

        if len(group_masks) == 0 or beadgroups[trial_bead] == len(group_masks):
            # Only DNA will be moving, create place-holder dummy coordinates
            p_coor_rot = np.zeros((0, 3))
        else:
            p_mask = group_masks[beadgroups[trial_bead]]
            p_ind_rot = mask2ind(p_mask)
            p_ind_fix = mask2ind(-(p_mask - 1))
            p_coor_rot = p_coor[p_ind_rot]
            p_coor_fix = p_coor[p_ind_fix]

        # generate a newly rotated model
        (d_coor[trial_bead:], xyz[:, trial_bead:], p_coor_rot
         ) = dna_move.beadRotate(d_coor[trial_bead - 1:], xyz[:, trial_bead - 1:],
                                 thetaXYZ, ARGS.n_soft, p_coor_rot)

        # store the rotated protein coordinates
        if beadgroups[trial_bead] < len(group_masks):
            p_coor[p_ind_rot] = p_coor_rot

        # verify the Rg_new < Rg_old * 1.01
        d_coor_old = np.copy(cg_dna.coor()[0])
        cg_dna.setCoor(np.array([(d_coor)]))  # update dna coordinates
        rg_new = cg_dna.calcrg(0)

        if rg_new < rg_old * 1.01:
            rg_pass = True
            print 'rg_old * 1.01 < rg_new: %f < %f' % (rg_old * 1.01, rg_new)
        else:
            rg_pass = False
            print 'rg_old * 1.01 > rg_new: %f > %f' % (rg_old * 1.01, rg_new)

        if rg_pass:
            # calculate the change in energy (dU) and the boltzman factor (p)
            (u, l) = dna_move.checkU(d_coor)
            Ub1 = dna_move.energyBend(lpl, u, l)

            # ~~~~ DNA interaction energy  ~~~~~~#
            (Uwca1, wca1) = dna_move.f_energy_wca(w, d_coor, wca0, trial_bead)

            U_T1 = Ub1 + Uwca1
            dU = U_T1 - U_T0

            with warnings.catch_warnings():
                warnings.filterwarnings('error')  # need this for np warnings
                try:
                    p = np.exp(-dU)
                except Warning:
                    if dU > 99:
                        p = 0
                        # s print 'energy was large, setting probability to 0'
                    elif dU < 0:
                        p = 1
                        # s print 'energy was negative, setting probability to
                        # 1'
                    else:
                        print 'Warning: ~~> unclear OverflowError <~~ dU = ', dU
                        print 'not sure where the error originated from'

            test = np.random.random()

            if p <= test:
                dna_pass = False
                # print 'step failed because of DNA energy'
            else:
                dna_pass = True

                # now check for collisions protein involved collisions
                if len(p_coor_rot) > 0:   # only if proteins were rotated
                    # ~~~~ Check for overlap, DNA-protein or protein-protein ~~~~~~#
                    d_coor_fix = d_coor[trial_bead:]
                    d_coor_rot = d_coor[:trial_bead]

                    # check for protein-protein overlap
                    if 1 == f_overlap2(p_coor_rot, p_coor_fix, pro_pro_test):
                        print 'Protein-Protein'
                        # print 'collision, set p=0'
                        collisionless = False

                    # print 'currently ignoring DNA-protein overlap'
                    # check for DNA-protein overlap
                    elif 1 == f_overlap2(p_coor_rot, d_coor_fix, dna_pro_test):
                        print 'Potein-DNA (rot-fix)'
                        # print 'collision, set p=0'
                        collisionless = False
                        print 'ignoring this for now'

                    elif 1 == f_overlap2(p_coor_fix, d_coor_rot, dna_pro_test):
                        print 'Potein-DNA (fix-rot)'
                        # print 'collision, set p=0'
                        collisionless = False
                    else:
                        collisionless = True

                    if not collisionless:
                        print 'failed because of collision'
                else:
                    collisionless = True  # no protein to collide with

        if rg_pass and dna_pass and collisionless:
            rg_old = rg_new
            n_from_reload += 1
            steps_from_0[n_accept] = n_from_reload + n_reload[-1]
            n_accept += 1                      # increment accept counter
            # cg_dna.setCoor(d_coor) # <-- DO NOT use setCoor, want uniuqe mem
            # cg_pro.setCoor(p_coor) # <-- DO NOT use setCoor, want uniuqe mem
            cg_pro.setCoor(np.array([(p_coor)]))  # update protein coordinates
            vecXYZ = np.copy(xyz)              # update dna orientations
            vecX_mol.setCoor(np.array([vecXYZ[0]]))  # independent of vecXYZ[0]
            vecY_mol.setCoor(np.array([vecXYZ[1]]))  # independent of vecXYZ[1]
            vecZ_mol.setCoor(np.array([vecXYZ[2]]))  # independent of vecXYZ[2]

            wca0 = np.copy(wca1)               # update DNA WCA energy
            U_T0 = U_T1                        # update total energy

            # print output regarding trial
            print "trial_bead(%3d) = %2d\t failed attempts = %2d" % (n_accept,
                                                                     trial_bead, fail_tally)
            fail_tally = 0                     # reset fail_tally

            # print out the Rg

            print cg_dna.calcrg(0)

            # write out the accepted configuration for go-back use
            if ARGS.goback > 0:
                # these are incremented by one because the original coordinates
                # are saved (that is not the case for aa_all)
                cg_dna.write_dcd_step(cg_dna_dcd_out, 0, n_written + 1)
                cg_pro.write_dcd_step(cg_pro_dcd_out, 0, n_written + 1)
                vecX_mol.write_dcd_step(vecX_dcd_out, 0, n_written + 1)
                vecY_mol.write_dcd_step(vecY_dcd_out, 0, n_written + 1)
                vecZ_mol.write_dcd_step(vecZ_dcd_out, 0, n_written + 1)

            # recover an all atom representation and save coordinates to a dcd
            # this requires re-inserting the aa-coordinates which takes added
            # time so only do when designated
            if 0 == n_accept % ARGS.n_dcd_write:
                # ~~recover aa-DNA~~
                error = dna_move.recover_aaDNA_model(cg_dna, aa_dna, vecXYZ,
                                                     all_beads, dna_bead_masks)
                # ~~recover aa-Protein~~
                dna_move.recover_aaPro_model(aa_pgroup_masks, cg_pgroup_masks,
                                             cg_pro, all_proteins, aa_pro)
                # ~~Combine aa Complete Structure~~
                aa_all.set_coor_using_mask(aa_pro, 0, aa_pro_mask)
                aa_all.set_coor_using_mask(aa_dna, 0, aa_dna_mask)
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written)

        else:
            # default ARGS.goback is -1 so this returns FALSE without user
            # input
            if fail_tally == ARGS.goback:
                i_goback = rewind(ARGS, n_accept, cg_dna_dcd_name,
                                  cg_dna, cg_pro_dcd_name, cg_pro, vecX_dcd_name,
                                  vecX_mol, vecY_mol, vecY_dcd_name, vecZ_mol,
                                  vecZ_dcd_name, vecXYZ)

                # revert dna coordinates
                cg_dna.setCoor(np.array([(d_coor_old)]))
                d_coor = np.copy(cg_dna.coor()[0])  # reset the dna coordinates

                # reset the reference energy
                (u, l) = checkU(d_coor)
                (Uwca0, wca0) = f_energy_wca(w, d_coor, wca0, 0)
                Ub0 = energyBend(lpl, u, l)
                U_T0 = Ub0 + Uwca0

                n_from_reload = 0
                n_reload.append(steps_from_0[i_goback - 1])
                fail_tally = 0  # reset the fail counter
            else:
                fail_tally += 1                 # increment bead reject counter
                # increment total reject counter
                n_reject += 1
                # revert dna coordinates
                cg_dna.setCoor(np.array([(d_coor_old)]))
                d_coor = np.copy(cg_dna.coor()[0])  # reset the dna coordinates

            p_coor = np.copy(cg_pro.coor()[0])  # reset the protein coordinates
            xyz = np.copy(vecXYZ)              # reset the dna orientations

            # save previous coordinates again
            if not ARGS.keep_unique:
                # ~~Write DCD step~~
                n_written += 1
                aa_all.write_dcd_step(aa_all_dcd_out, 0, n_written)
                cg_dna.write_dcd_step(cg_dna_all_dcd_out, 0, n_written + 1)
                cg_pro.write_dcd_step(cg_pro_all_dcd_out, 0, n_written + 1)

    aa_all.close_dcd_write(aa_all_dcd_out)

    cg_dna.close_dcd_write(cg_dna_dcd_out)  # uncomment if wanting to keep
    # os.remove(timestr + 'cg_dna.pdb') # remove/comment to keep the cg dna coor
    # os.remove(cg_dna_dcd_name) # remove/comment to keep the cg dna coor
    # cg_pro.close_dcd_write(cg_pro_dcd_out) #uncomment if wanting to keep
    os.remove(timestr + 'cg_pro.pdb')  # remove/comment to keep the cg pro coor
    os.remove(cg_pro_dcd_name)  # remove/comment to keep the cg pro coor
    os.remove(vecX_dcd_name)
    os.remove(vecY_dcd_name)
    os.remove(vecZ_dcd_name)

    if ARGS.goback > 0:
        np.savetxt(timestr + 'n_from_0.txt', steps_from_0, fmt='%d')

    print "accepted %d moves" % n_accept
    print "rejected %d moves" % n_reject

    # print n_reload
    # print steps_from_0


def main():

    # set the DNA properties parameters:
    lp = 530.     # persistence length  (lp = 530A)
    # defined in dna_mc
    # w = 46        # effective measure of dsDNA chain width in A (w = 46A)
    # DOI 10.1021/ma201277e used w=46, lp=530

    print 'loading pdb'
    all_atom_pdb = ARGS.pdb
    dna_resids = []
    pro_groups = []

    if ARGS.pdb == 'new_c11_tetramer.pdb':
        '''The C11 nucleosome tetramer with all the protein tails'''
        dna_segnames = ['DNA1', 'DNA2']
        dna_resids.append([1, 693])
        dna_resids.append([693, 1])
        # continuous flexible residues on the first DNA strand
        # recall that range(a, b) excludes upper lim: [a, b)
        flex_resids = [range(1, 31), range(167, 198), range(334, 365),
                       range(501, 532), range(667, 694)]
        pro_groups.append(['A0', 'B0', 'C0', 'D0',
                           'E0', 'F0', 'G0', 'H0'])
        pro_groups.append(['A1', 'B1', 'C1', 'D1',
                           'E1', 'F1', 'G1', 'H1'])
        pro_groups.append(['M1', 'N1', 'O1', 'P1',
                           'Q1', 'R1', 'S1', 'T1'])
        pro_groups.append(['M0', 'N0', 'O0', 'P0',
                           'Q0', 'R0', 'S0', 'T0'])
        ARGS.theta_max = [10, 10, 10, 10, 10]

    elif ARGS.pdb == 'new_c11h5.pdb':
        '''The C11 nucleosome tetramer with the gH5 linker'''
        dna_segnames = ['DNA1', 'DNA2']
        dna_resids.append([1, 693])
        dna_resids.append([693, 1])
        # continuous flexible residues on the first DNA strand
        # s flex_resids = [range(1, 31), range(167, 198), range(334, 365),
        # s     range(501, 532), range(667, 694)]
        small_range = [range(1, 10), range(177, 180), range(344, 352),
                       range(513, 517), range(680, 694)]
        inter_range = [range(1, 20), range(172, 189), range(339, 357),
                       range(507, 524), range(674, 694)]
        large_range = [range(1, 31), range(167, 198), range(334, 365),
                       range(501, 532), range(667, 694)]
        xlarge_range = [range(1, 36), range(162, 203), range(329, 370),
                        range(496, 537), range(662, 699)]

        if ARGS.size == 'large':
            flex_resids = large_range
        elif ARGS.size == 'med':
            flex_resids = inter_range
        elif ARGS.size == 'small':
            flex_resids = small_range
        elif ARGS.size == 'x-lrg':
            flex_resids = xlarge_range
        else:
            print 'unknown selection for size: ', size

        if ARGS.regions == '1-5':         # 1-5 option
            pass                          # flex_resids = flex_resids
        elif ARGS.regions == '234':       # 2-4 option
            flex_resids[0] = flex_resids[4] = []
        elif ARGS.regions == '2_4':       # 2, 4 option
            flex_resids[0] = flex_resids[4] = flex_resids[2] = []
        elif ARGS.regions == '23_':       # 2-3 option
            flex_resids[0] = flex_resids[4] = flex_resids[3] = []
        elif ARGS.regions == '2__':       # 2 option
            flex_resids[0] = flex_resids[
                4] = flex_resids[3] = flex_resids[2] = []
        elif ARGS.regions == '_3_':       # 3 option
            flex_resids[0] = flex_resids[
                4] = flex_resids[3] = flex_resids[1] = []
        else:
            print 'unknown selection for flex: ', ARGS.regions

        # manually selecting flex-resids
        # s flex_resids = [range(18,33), range(182,199), [], range(496,508),
        # range(668,677)]
        flex_resids = [[], [], [], [], range(668, 677)]

        print ('using regions ' + ARGS.regions + ' with a ' + ARGS.size +
               ' size -> flex_resids = \n', flex_resids)

        pro_groups.append(['A0', 'B0', 'C0', 'D0',
                           'E0', 'F0', 'G0', 'H0', 'H5N1'])
        pro_groups.append(['A1', 'B1', 'C1', 'D1',
                           'E1', 'F1', 'G1', 'H1', 'H5N2'])
        pro_groups.append(['M1', 'N1', 'O1', 'P1',
                           'Q1', 'R1', 'S1', 'T1', 'H5N3'])
        pro_groups.append(['M0', 'N0', 'O0', 'P0',
                           'Q0', 'R0', 'S0', 'T0', 'H5N4'])
        ARGS.theta_max = [2, 5, 2, 5, 2]

    elif ARGS.pdb == 'new_dsDNA60.pdb':
        # linker dna file
        dna_segnames = ['DNA1', 'DNA2']
        dna_resids.append([1, 60])  # DNA base pairing
        dna_resids.append([120, 61])  # DNA base pairing
        flex_resids = [range(16, 45)]
        ARGS.theta_max = [25, 25, 25, 25, 25]
    elif ARGS.pdb == 'new_dsDNA.pdb':
        # linker dna file
        dna_segnames = ['DNA1', 'DNA2']
        dna_resids.append([1, 30])  # DNA base pairing
        dna_resids.append([30, 1])  # DNA base pairing
        flex_resids = [range(1, 31)]
    else:
        print "\n~~~ ERROR, unknow pdb file input ~~~\n"

    tic = time.time()

    (cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, trialbeads, beadgroups,
     move_masks, all_beads, dna_bead_masks, aa_pgroup_masks, cg_pgroup_masks,
     all_proteins, aa_all, aa_pro_mask, aa_dna_mask, bp_per_bead
     ) = dna_move.get_cg_parameters(
        ARGS, flex_resids, pro_groups, dna_resids, dna_segnames)

    toc = time.time() - tic
    print 'Total coarse-grain time = %f seconds' % toc

    tic = time.time()
    driven_dna_mc(ARGS, cg_dna, aa_dna, cg_pro, aa_pro, vecXYZ, lp, trialbeads,
                  beadgroups, move_masks, all_beads, dna_bead_masks, aa_pgroup_masks,
                  cg_pgroup_masks, all_proteins, aa_all, aa_pro_mask, aa_dna_mask)
    toc = time.time() - tic
    print 'run time =', toc, 'seconds'

    if None == ARGS.Llp and ARGS.flex:
        write_flex_resids(all_beads, trialbeads, ARGS)

if __name__ == "__main__":

    import argparse

    ''' wanted to implement this for error handling but never did'''
    # import logging
    # LOGGER = logging.getLogger(__name__) #add module name manually
    # import sys
    # if '-v' in sys.argv:
    #     logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
    # else:
    #     logging.basicConfig()

    ARGS = dna_move.parse()  # this makes ARGS global

    main()
    print '\nFinished %d successful DNA moves! \n\m/ >.< \m/' % ARGS.nsteps
