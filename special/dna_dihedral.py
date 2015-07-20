#!/usr/bin/python
#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: methods for calculating DNA dihedral angles and generating plots
# Created: 22 September 2014
#
# $Id$
#
# 0000000011111111112222222222333333333344444444445555555555666666666677777777778
# 2345678901234567890123456789012345678901234567890123456789012345678901234567890

import sassie.sasmol.sasmol as sasmol
import numpy
import matplotlib
import matplotlib.pyplot as plt
import logging
import sys

import sassie.simulate.monte_carlo.monomer.dihedral_monte_carlo as dmc
#import dihedral_monte_carlo as dmc
import sassie.simulate.monte_carlo.monomer.dihedral_rotate as dr
#import dihedral_rotate as dr


def plot_me(x_angle, y_angle, x_angle_label, y_angle_label):

    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    ax1 = plt.subplot(111)
    ax1.scatter(x_angle, y_angle)

    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-180, 180])

    ax1.set_xlabel(x_angle_label)
    ax1.set_ylabel(y_angle_label)
    plt.title(y_angle_label + ' vs ' + x_angle_label)

    plt.show()


def plot_dna_dihedral_grid(all_angles):

    plt.ion()
    (n_frames, n_bps, n_angles) = all_angles.shape

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}

    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'

        zero_to_threesixty(all_angles[frame])
        plt.plot(all_angles[frame][0:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off', labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        plt.plot(all_angles[frame][0:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off', labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        plt.plot(all_angles[frame][:, angles[x]],
                 all_angles[frame][:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        plt.plot(all_angles[frame][:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        lines = []
        labels = []
        plt.plot(all_angles[frame][:-1, angles[x]],
                 all_angles[frame][1:, angles[y]], 'o')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')

        plt.suptitle('Starting Structure')
        plt.savefig('DNA-dihedrals_0steps.pdf', bbox_ingches='tight')
        plt.savefig('DNA-dihedrals_0steps.png', bbox_ingches='tight')

        plt.draw()

    plt.show()


def plot_dna_compare_dihedral(n0_all_angles, nm_all_angles, m__all_angles):
    angles = {'alpha': 0, 'beta': 1, 'gamma': 2, 'delta': 3,
              'epsilon': 4, 'zeta': 5, 'chi': 6}
    zero_to_threesixty(n0_all_angles)
    zero_to_threesixty(nm_all_angles)
    zero_to_threesixty(m_all_angles)

    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    plt.title('Scatter plots of selected torsional angles.')

    ax1 = plt.subplot(331)
    x = 'zeta'
    y = 'alpha'
    ax1.scatter(n0_all_angles[angles[x]][
                0:-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax1.scatter(nm_all_angles[angles[x]][
                0:-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax1.scatter(m__all_angles[angles[x]][
                0:-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y + ' + 1')
    grid_format(ax1)

    ax2 = plt.subplot(332)
    x = 'zeta'
    y = 'beta'
    ax2.scatter(n0_all_angles[angles[x]][
                0:-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax2.scatter(nm_all_angles[angles[x]][
                0:-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax2.scatter(m__all_angles[angles[x]][
                0:-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax2.set_xlabel(x)
    ax2.set_ylabel(y + ' + 1')
    grid_format(ax2)

    all_angles = n0_all_angles

    ax3 = plt.subplot(333)
    x = 'zeta'
    y = 'epsilon'
    ax3.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax3.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax3.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax3.set_xlabel(x)
    ax3.set_ylabel(y)
    grid_format(ax3)

    ax4 = plt.subplot(334)
    x = 'alpha'
    y = 'gamma'
    ax4.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax4.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax4.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax4.set_xlabel(x)
    ax4.set_ylabel(y)
    grid_format(ax4)

    ax5 = plt.subplot(335)
    x = 'zeta'
    y = 'chi'
    ax5.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax5.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax5.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax5.set_xlabel(x)
    ax5.set_ylabel(y)
    grid_format(ax5)

    ax6 = plt.subplot(336)
    x = 'delta'
    y = 'chi'
    ax6.scatter(n0_all_angles[angles[x]], n0_all_angles[
                angles[y]], c='b', marker='x', s=80)
    ax6.scatter(nm_all_angles[angles[x]], nm_all_angles[
                angles[y]], c='r', marker='o', s=25)
    ax6.scatter(m__all_angles[angles[x]], m__all_angles[
                angles[y]], c='g', marker='s', s=25)
    ax6.set_xlabel(x)
    ax6.set_ylabel(y)
    grid_format(ax6)

    ax7 = plt.subplot(337)
    x = 'zeta'
    y = 'zeta'
    ax7.scatter(n0_all_angles[angles[x]][
                :-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    ax7.scatter(nm_all_angles[angles[x]][
                :-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    ax7.scatter(m__all_angles[angles[x]][
                :-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    ax7.set_xlabel(x)
    ax7.set_ylabel(y + ' + 1')
    grid_format(ax7)

    ax8 = plt.subplot(338)
    x = 'epsilon'
    y = 'epsilon'
    ax8.scatter(n0_all_angles[angles[x]][
                :-1], n0_all_angles[angles[y]][1:], c='b', marker='x', s=80, label="starting")
    ax8.scatter(nm_all_angles[angles[x]][
                :-1], nm_all_angles[angles[y]][1:], c='r', marker='o', s=25, label="final")
    ax8.scatter(m__all_angles[angles[x]][
                :-1], m__all_angles[angles[y]][1:], c='g', marker='s', s=25, label="minimized")
    ax8.set_xlabel(x)
    ax8.set_ylabel(y + ' + 1')
    grid_format(ax8)
    ax8.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    # ax9 = plt.subplot(339)
    # x = 'zeta'
    # y = 'P'
    # ax9.scatter(n0_all_angles[angles[x]][:-1],n0_all_angles[angles[y]][1:], c='b', marker='x', s=80)
    # ax9.scatter(nm_all_angles[angles[x]][:-1],nm_all_angles[angles[y]][1:], c='r', marker='o', s=25)
    # ax9.scatter(m__all_angles[angles[x]][:-1],m__all_angles[angles[y]][1:], c='g', marker='s', s=25)
    # ax9.set_xlim([0, 360])
    # ax9.set_ylim([0, 360])
    # ax9.set_xlabel(x)
    # ax9.set_ylabel(y)

    plt.show()


def plot3d_dihedral(min_angles, raw_angles):
    from mayavi import mlab

    logging.debug('testing')
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape
    (n_thetas_r, n_frames_r, n_bps_r, n_angles_r) = raw_angles.shape


def plot_dna_min_dihedral(min_angles, raw_angles):
    logging.debug('testing')
    plt.ion()
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape
    (n_thetas_r, n_frames_r, n_bps_r, n_angles_r) = raw_angles.shape

    # assert that raw_angle parameters match the min parameters

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    step = 1
    frame_0 = step - 1

    for frame in xrange(frame_0, n_frames, step):
        # for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'
        for tm in xrange(n_thetas):
            # only need to do this once
            zero_to_threesixty(min_angles[tm, frame])
            zero_to_threesixty(raw_angles[tm, frame])
            plt.plot(raw_angles[tm, frame][0:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off', labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][0:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off', labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:, angles[x]],
                     raw_angles[tm, frame][:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        lines = []
        labels = []
        for tm in xrange(n_thetas):
            plt.plot(raw_angles[tm, frame][:-1, angles[x]],
                     raw_angles[tm, frame][1:, angles[y]], 'r+')
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')
        plt.legend(['raw', 'min'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        # plt.subplot(339)
        # x = r'$\zeta$'
        # y = r'$\delta$'
        # lines = []
        # labels = []
        # for tm in xrange(n_thetas):
        # plt.plot(all_angles[tm,frame][:, angles[x]], all_angles[tm,frame][:, angles[y]], 'gx')
        # labels.append('run %d' % theta_max[tm])
        # # plt.tick_params(labelleft='off',labelright='on')
        # # plt.xlabel(x)
        # # plt.ylabel(y)
        # # plt_format()
        # plt.tick_params(labelleft='off',labelbottom='off')
        # plt.xlim([-999, 999])
        # plt.ylim([-999, 999])
        # plt.legend(labels)

        # plt.suptitle('Starting Structure')
        # plt.savefig('DNA-dihedrals_0steps.pdf', bbox_ingches='tight')

        plt.suptitle(
            '%d Steps: scatter plots of selected torsional angles' % ((frame + 1) * 100))
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.pdf' %
                    ((frame + 1) * 100), bbox_inches='tight')
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.png' %
                    ((frame + 1) * 100), bbox_inches='tight')
        plt.draw()

    plt.show()


def plot_dna_dihedral(min_angles):
    logging.debug('testing')
    plt.ion()
    (n_thetas, n_frames, n_bps, n_angles) = min_angles.shape

    # assert that raw_angle parameters match the min parameters

    angles = {r'$\alpha$': 1, r'$\beta$': 2, r'$\gamma$': 3, r'$\delta$': 4,
              r'$\epsilon$': 5, r'$\zeta$': 6, r'$\chi$': 7}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    step = 1
    frame_0 = step - 1

    for frame in xrange(frame_0, n_frames, step):
        # for frame in xrange(n_frames):
        plt.figure()

        plt.subplot(331)
        x = r'$\zeta$'
        y = r'$\alpha$'
        for tm in xrange(n_thetas):
            # only need to do this once
            zero_to_threesixty(min_angles[tm, frame])
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')

        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelbottom='off', labeltop='on')

        plt.subplot(332)
        x = r'$\zeta$'
        y = r'$\beta$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][0:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labeltop='on')

        plt.subplot(333)
        x = r'$\zeta$'
        y = r'$\epsilon$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(
            labelleft='off', labelbottom='off', labeltop='on', labelright='on')

        plt.subplot(334)
        x = r'$\alpha$'
        y = r'$\gamma$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelbottom='off')

        plt.subplot(335)
        x = r'$\zeta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off')

        plt.subplot(336)
        x = r'$\delta$'
        y = r'$\chi$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:, angles[x]],
                     min_angles[tm, frame][:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y)
        plt_format()
        plt.tick_params(labelleft='off', labelbottom='off', labelright='on')

        plt.subplot(337)
        x = r'$\zeta$'
        y = r'$\zeta$'
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()

        plt.subplot(338)
        x = r'$\epsilon$'
        y = r'$\epsilon$'
        lines = []
        labels = []
        for tm in xrange(n_thetas):
            plt.plot(min_angles[tm, frame][:-1, angles[x]],
                     min_angles[tm, frame][1:, angles[y]], 'gx')
        plt.xlabel(x)
        plt.ylabel(y + ' + 1')
        plt_format()
        plt.tick_params(labelleft='off')
        plt.legend(['min', 'raw'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        # plt.subplot(339)
        # x = r'$\zeta$'
        # y = r'$\delta$'
        # lines = []
        # labels = []
        # for tm in xrange(n_thetas):
        # plt.plot(all_angles[tm,frame][:, angles[x]], all_angles[tm,frame][:, angles[y]], 'gx')
        # labels.append('run %d' % theta_max[tm])
        # # plt.tick_params(labelleft='off',labelright='on')
        # # plt.xlabel(x)
        # # plt.ylabel(y)
        # # plt_format()
        # plt.tick_params(labelleft='off',labelbottom='off')
        # plt.xlim([-999, 999])
        # plt.ylim([-999, 999])
        # plt.legend(labels)

        # plt.suptitle('Starting Structure')
        # plt.savefig('DNA-dihedrals_0steps.pdf', bbox_ingches='tight')

        plt.suptitle(
            '%d Steps: scatter plots of selected torsional angles' % ((frame + 1) * 10))
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.pdf' %
                    ((frame + 1) * 10), bbox_inches='tight')
        plt.savefig('DNA-dihedrals_best_%dsteps-tz.png' %
                    ((frame + 1) * 10), bbox_inches='tight')
        plt.draw()

    plt.show()


def plt_format():

    plt.xlim([0, 360])
    plt.ylim([0, 360])
    plt.xticks([0, 60, 120, 180, 240, 300, 360])
    plt.yticks([0, 60, 120, 180, 240, 300, 360])


def grid_format(ax):

    ax.set_xlim([0, 360])
    ax.set_ylim([0, 360])
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_yticks([0, 60, 120, 180, 240, 300, 360])


def zero_to_threesixty(all_angles):
    for (i, angle_list) in enumerate(all_angles):
        for (j, angle) in enumerate(angle_list):
            if angle < 0:
                all_angles[i][j] = angle + 360


def load_angles(file_name):
    in_file = open(file_name, 'r')
    i_frame = 0
    all_frames = []
    for line in in_file.readlines():
        if line.find("#", 0) == 0:
            if line.find("frame") > 0:
                # start a new frame
                i_frame += 1
                if i_frame > 1:
                    all_frames.append(frame)
                frame = []
                continue
            continue
        frame.append(numpy.fromstring(line[2:], dtype='float64', sep='\t'))

    all_frames.append(frame)
    return numpy.array(all_frames)


def get_angles(dcd_file_name, pdb_file_name, first_last_resid, flex_file, nf=0):
    from dna.ds_dna_monte_carlo import read_flex_resids

    txtOutput = []

    out_file = open(dcd_file_name[:-4] + '.ddat', 'w')

    # read flex file
    # flex_file = pdb_file_name[:-3] + 'flex'
    (segnames, flex_resids) = read_flex_resids(flex_file)

    segname = segnames[0]

    numranges = 0
    resid_old = -999.999  # any number that is not an int
    for resid_new in flex_resids[:, 0]:
        if resid_new != resid_old + 1:
            numranges += 1
        resid_old = resid_new

    res_low = [flex_resids[0, 0]]
    n_cont = [len(flex_resids[:, 0])]
    flexible_residues = dmc.get_flexible_residues(numranges, res_low, n_cont)

    molecule_type = "dna"

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file_name)

    residue_rotation_indices, residue_rotation_mask = dmc.get_rotation_indices(
        mol, molecule_type, flexible_residues, txtOutput)

    # print "residue_rotation_indices = ",residue_rotation_indices

    # first_last_resid = [flexible_residues[0], flexible_residues[-1]]

    alpha = []
    beta = []
    gamma = []
    delta = []     # <--- significant
    epsilon = []
    zeta = []
    chi = []       # <--- significant

    all_alpha = []
    all_beta = []
    all_gamma = []
    all_delta = []     # <--- significant
    all_epsilon = []
    all_zeta = []
    all_chi = []       # <--- significant

    dcdfile = mol.open_dcd_read(dcd_file_name)
    if nf == 0:
        nf = dcdfile[2]  # number of frames

    # reslist = range(first_last_resid[0],first_last_resid[1]+1)
    # print 'reslist = ',reslist

    st1 = "# base\talpha\tbeta\tgamma\tdelta\tepsilon\tzeta\tchi\n"
    out_file.write("%s" % (st1))

    res_mol = sasmol.SasMol(0)

    for i in xrange(nf):

        mol.read_dcd_step(dcdfile, i)
        coor = mol.coor()
        out_file.write("# frame %d\n" % (i + 1))

        for res in xrange(n_cont[0]):
            q0 = res_low[0] + res
            error, res_mask = mol.get_subset_mask(
                "(resid[i] == " + str(q0) + ")")
            mol.copy_molecule_using_mask(res_mol, res_mask, 0)
            st = '%s %d' % (res_mol.resname()[0][0], q0)

            # get the indices and mask for this residue
            indices = residue_rotation_indices[q0]
            this_mask = numpy.array(residue_rotation_mask[q0])

            # alpha
            alpha_angle = dr.measure(coor, indices, "alpha", this_mask, q0,
                                     first_last_resid, molecule_type)
            str_alpha_angle = '%.1f' % alpha_angle
            all_alpha.append(alpha_angle)
            st = st + '\t' + str_alpha_angle
            logging.debug('alpha: ' + str_alpha_angle)

            # beta
            beta_angle = dr.measure(coor, indices, "beta", this_mask, q0,
                                    first_last_resid, molecule_type)
            str_beta_angle = '%.1f' % (beta_angle)
            all_beta.append(beta_angle)
            st = st + '\t' + str_beta_angle
            logging.debug('beta: ' + str_beta_angle)

            # gamma
            gamma_angle = dr.measure(coor, indices, "gamma", this_mask, q0,
                                     first_last_resid, molecule_type)
            str_gamma_angle = '%.1f' % (gamma_angle)
            all_gamma.append(gamma_angle)
            st = st + '\t' + str_gamma_angle
            logging.debug('gamma: ' + str_gamma_angle)

            # delta
            delta_angle = dr.measure(coor, indices, "delta", this_mask, q0,
                                     first_last_resid, molecule_type)
            str_delta_angle = '%.1f' % (delta_angle)
            all_delta.append(delta_angle)
            st = st + '\t' + str_delta_angle
            logging.debug('delta: ' + str_delta_angle)

            # epsilon
            epsilon_angle = dr.measure(coor, indices, "epsilon", this_mask, q0,
                                       first_last_resid, molecule_type)
            str_epsilon_angle = '%.1f' % (epsilon_angle)
            all_epsilon.append(epsilon_angle)
            st = st + '\t' + str_epsilon_angle
            logging.debug('epsilon: ' + str_epsilon_angle)

            # zeta
            zeta_angle = dr.measure(coor, indices, "zeta", this_mask, q0,
                                    first_last_resid, molecule_type)
            str_zeta_angle = '%.1f' % (zeta_angle)
            all_zeta.append(zeta_angle)
            st = st + '\t' + str_zeta_angle
            logging.debug('zeta: ' + str_zeta_angle)

            # chi
            chi_angle = dr.measure(coor, indices, "chi", this_mask, q0,
                                   first_last_resid, molecule_type)
            str_chi_angle = '%.1f' % (chi_angle)
            all_chi.append(chi_angle)
            st = st + '\t' + str_chi_angle
            logging.debug('alpha: ' + str_alpha_angle)

            out_file.write("%s\n" % (st))

    all_angles = [all_alpha, all_beta, all_gamma, all_delta, all_epsilon,
                  all_zeta, all_chi]

    out_file.close()

    print '\nsuccessfully finished calculating dihedrals \_/'
    return all_angles


def main():
    NotImplemented

if __name__ == '__main__':

    if '-v' in sys.argv:
        logging.basicConfig(level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

        main()
