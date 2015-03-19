#!/usr/bin/env python
# Auther: Steven C. Howell
# Purpose: Run crysol in parrallel with consolidated output
# Created: 10/09/2014
# $Id$

import sys
import os
import errno
import os.path as op
import subprocess
import logging
import cmd
import shutil
import time
import glob
# import sassie_1_na.drivers.my_crysol_driver as crysol

LOGGER = logging.getLogger(__name__) #add module name manually


class MainError(Exception):
    pass

def parse():
    ''' Returns arguments in parser'''
    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'Generate modified DNA or DNA-protein structures'
        #epilog = 'no epilog found'
        )
    parser.add_argument('-n', '--ncpu', default=4, type=int,
                        help = 'number of cpus to use for calculation'
                        )
    parser.add_argument('-r', '--runname', default='run0', type=str,
                        help = 'folder to put output files'
                        )
    parser.add_argument('-p', '--pdb', default='new_c11h5.pdb', type=str,
                        help = 'pdb template for structures to be calculated'
                        )
    parser.add_argument('-d', '--dcd', default='new_c11h5.dcd', type=str,
                        help = 'dcd file containing structures to calculate'
                        )
    parser.add_argument('-s', '--sleep', default=60, type=int,
                        help = 'time between checks for crysol to finish'
                        )
    parser.add_argument('-dv', '--driver', default='/home/schowell/data/code/pylib/sassie_1_na/drivers/my_crysol_driver.py', 
                        type=str, help = 'path to crysol driver'
                        )    
    parser.add_argument('-lm', '--maxh', default=4, type=int,
                        help = 'maximum order of harmonics for crysol'
                        )
    parser.add_argument('-fb', '--fib', default=4, type=int,
                        help = 'order of Fibonacci grid for crysol'
                        )
    parser.add_argument('-ns', '--numpoints', default=4, type=int,
                        help = 'number of points in crysol output'
                        )
    parser.add_argument('-sm', '--maxs', default=4, type=float,
                        help = 'maxs s (or q) of crysol output'
                        )
    return parser.parse_args()

def mkdir_p(path):
    '''
    make directory recursively
    adapted from http://stackoverflow.com/questions/600268/
    '''
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def append_bk(folder):
    new_folder = folder + '_BK'
    if os.path.exists(new_folder):
        append_bk(new_folder)
    else:
        shutil.move(folder, new_folder)
        print 'moved %s to %s' % (folder, new_folder)

class folder_exists(cmd.Cmd):
    
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(0/1/2)> '
        
    def do_move(self, arg):
        append_bk(self.runname)
        return True
    
    def help_move(self):
        print '-- move run folder to run_BK'

    def do_replace(self, arg):
        print 'removing run folder'
        shutil.rmtree(self.runname)
        return True
    
    def help_replace(self):
        print 'remove and replace run folder'
    
    def do_quit(self, arg):
        print 'exiting program'
        sys.exit(1)
        
    def help_quit(self):
        print '-- terminates the application'

    def default(self, arg):
        print 'invalid selection, please select: 0/1/2'
            

    #shortcuts
    do_0 = do_quit
    do_1 = do_move
    do_2 = do_replace
    help_0 = help_quit
    help_1 = help_move
    help_2 = help_replace

class cd:
    """
    Context manager for changing the current working directory
    http://stackoverflow.com/questions/431684
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def tail(f, n=10):
    '''
    return the last n lines of f
    adapted from: http://stackoverflow.com/questions/136168
    '''
    tail_str = 'tail -n %s %s' % (str(n), f)
    stdin,stdout = os.popen2(tail_str)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines[:]

def collect_crysol(sub_dirs, runname, sleep):
    out_dir = os.getcwd() + '/' + runname + '/crysol'
    mkdir_p(out_dir)
    
    n_out_files = 1
    
    for (i, sub_dir) in enumerate(sub_dirs):
        logging.debug('waiting for %s' % sub_dir)
        with cd(sub_dir):
            # read the end of the file test.out
            # out_file = 'par_crysol_%02d.out' % (i+1)
            out_files = glob.glob('*out')
            not_done = True
            while not_done:
                time.sleep(sleep)
                for (j, out_file) in enumerate(out_files):
                    if tail(out_file, 1) == ['GCRYSOL IS DONE\n']:
                        not_done = False

            print '%s is complete' % sub_dir
            sub_files = glob.glob('*/crysol/*int')
            more_sub_files = glob.glob('crysol/*int')
            n_sub_files = len(sub_files)
            n_more_sub_files = len(more_sub_files)
            if n_sub_files > 0 and n_more_sub_files > 0:
                logging.warning(('found crysol output in both "./*/crysol/" '
                    '(%d files) and "./crysol/" (%d files), order may be '
                    'unexpected') % (n_sub_files, n_more_sub_files) )
            for another_file in more_sub_files:
                sub_files.append(another_file)
                
            error = sub_files.sort()
            for (j, sub_file) in enumerate(sub_files):
                logging.debug('moving %s' % sub_file)
                file_name = sub_file[:-4]
                new_name = runname + '_' + str(n_out_files + j).zfill(5)
                logging.debug('moving %s%s.int to %s/%s.int' % (sub_dir, 
                                    file_name, out_dir, new_name))
                os.system('mv %s.int %s/%s.int' % (file_name, out_dir, new_name))
                os.system('mv %s.log %s/%s.log' % (file_name, out_dir, new_name))                
                
            # os.system('mv *out *dcd *pdb ../')
            os.system('mv *out ../')
            n_out_files += len(sub_files)
            
        shutil.rmtree(sub_dir)

def iterate_crysol(inputs, sub_dirs, dcd_file_names):
    if os.path.exists('/share/apps/bin/python'):
        #gibbs
        python = '/share/apps/bin/python'
    elif os.path.exists('/share/apps/local/bin/python/bin/python'):
        #entropy
        python = '/share/apps/local/bin/python/bin/python'
    elif os.path.exists('/usr/bin/python'):
        #my machines
        python = '/usr/bin/python'
    else:
        #others
        python = '/usr/bin/env python'
    
    # sub_inputs = crysol.inputs()
    # sub_inputs.runname   = 'par'
    # sub_inputs.pdbpath   = './'
    # sub_inputs.dcdpath   = './'
    # sub_inputs.pdbfile   = inputs.pdb
    # sub_inputs.maxh      = inputs.maxh
    # sub_inputs.fib       = inputs.fib
    # sub_inputs.maxs      = inputs.maxs
    # sub_inputs.numpoints = inputs.numpoints
    _, driver = os.path.split(inputs.driver)
    all_run_str = '%s %s -r par -pp ./ -p %s -dp ./ -lm %d -ns %d -sm %f -fb %d' % (python, driver, inputs.pdb, inputs.maxh, inputs.numpoints , inputs.maxs, inputs.fib)
    for (i, sub_dir) in enumerate(sub_dirs):
        sub_run_str = ' -d %s > par_crysol_%02d.out  &' % (dcd_file_names[i], i+1)
        run_str = all_run_str + sub_run_str
        with cd(sub_dir):
            os.system('cp %s ./' % inputs.driver)
            os.system(run_str)

def split_dcd(inputs):
    import sassie.sasmol.sasmol as sasmol
    import numpy as np
    inputs.out_dir = inputs.runname + '/crysol'
    if os.path.exists(inputs.out_dir):
        print 'WARNING: run folder exists (%s), moving it\n' % inputs.out_dir
        append_bk(inputs.out_dir)
        # print 'select one of the following (0/1/2): quit / move / replace'
        # folder = folder_exists()
        # folder.runname = inputs.out_dir
        # result = folder.cmdloop()
    else:
        print 'created new run folder: %s' % inputs.out_dir
    mkdir_p(inputs.out_dir)
    
    mol = sasmol.SasMol(0)
    mol.read_pdb(inputs.pdb)
    
    # mol.read_dcd(inputs.dcd)
    dcd_file = mol.open_dcd_read(inputs.dcd)
    total_frames = dcd_file[2]
    n_atoms = dcd_file[1]
    copy_mask = np.ones(n_atoms, dtype=np.int32)
    
    if inputs.ncpu < 0:
        print 'ncpu: %d < 0,   using |%d| = %d instead' % (inputs.ncpu, 
                                                inputs.ncpu, abs(inputs.ncpu) )
        inputs.ncpu = abs(inputs.ncpu)
    n_frames_sub = total_frames/inputs.ncpu
    last_frame = 0
    sub_dirs = []
    dcd_file_names = []
    for cpu in xrange(1, inputs.ncpu+1):
        sub_dir = inputs.out_dir + '/sub' + str(cpu).zfill(2) + '/'
        sub_dirs.append(sub_dir)
        mkdir_p(sub_dir)
        os.system('cp %s %s' % (inputs.pdb, sub_dir))
        sub_mol = sasmol.SasMol(0)
        mol.copy_molecule_using_mask(sub_mol, copy_mask, 0)
        with cd(sub_dir):
            if cpu == inputs.ncpu:
                n_frames_sub = n_frames_sub + total_frames % inputs.ncpu
            dcd_out_name = 'sub' + str(cpu).zfill(2) + '.dcd'
            dcd_file_names.append(dcd_out_name)
            first = last_frame
            last = last_frame + n_frames_sub
            dcd_out_file = sub_mol.open_dcd_write(dcd_out_name)
            for (i, frame) in enumerate(xrange(first, last)):
                sub_mol.read_dcd_step(dcd_file, frame)
                sub_mol.write_dcd_step(dcd_out_file, 0, i+1)
                
            sub_mol.close_dcd_write(dcd_out_file)

        del sub_mol
            
        last_frame += n_frames_sub
    print 
    return sub_dirs, dcd_file_names

def main(inputs):

    if inputs.debug:
        logging.basicConfig(filename='%s.log' % __file__[:-3], level=logging.DEBUG)
    else:
        logging.basicConfig()

    #check the input
    assert os.path.exists(inputs.pdb), 'ERROR: "%s" does not exist' % inputs.pdb 
    assert os.path.exists(inputs.dcd), 'ERROR: "%s" does not exist' % inputs.dcd
    assert os.path.exists(inputs.driver), 'ERROR: "%s" does not exist' % inputs.driver 

    #break dcd into N dcds with a folder for each
    sub_dirs, dcd_file_names = split_dcd(inputs)
    
    #run crysol instance on each folder
    iterate_crysol(inputs, sub_dirs, dcd_file_names)
    
    #collect the results
    collect_crysol(sub_dirs, inputs.runname, inputs.sleep)

    print 'finished all crysol calculations\n     \m/ >.< \m/ '

if __name__ == '__main__':
    import argparse
    
    if '-v' in sys.argv:
        logging.basicConfig(filename='%s.log' %__file__[:-3], level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    ARGS = parse()    
    main(ARGS)
