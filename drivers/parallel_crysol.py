#!/usr/bin/env python
# Auther: Steven C. Howell
# Purpose: Run crysol in parrallel with consolidated output
# Created: 10/09/2014
# $Id: parallel_crysol.py,v 1.2 2014-10-09 23:28:44 schowell Exp $

import sys
import os
import os.path as op
import subprocess
import logging
import cmd
import shutil
import time
import glob

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
    return parser.parse_args()

def append_bk(folder):
    new_folder = folder + '_BK'
    if os.path.exists(new_folder):
        append_bk(new_folder)
    else:
        shutil.move(ARGS.runname,new_folder)
        print 'moved %s to %s' % (ARGS.runname,new_folder)

class folder_exists(cmd.Cmd):
    
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(0/1/2)> '
        
    def do_move(self, arg):
        append_bk(ARGS.runname)
        return True
    
    def help_move(self):
        print '-- move run folder to run_BK'

    def do_replace(self, arg):
        print 'removing run folder'
        shutil.rmtree(ARGS.runname)
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

def collect_crysol(sub_dirs):
    outdir = ARGS.runname + '/crysol'
    os.mkdir(outdir)
    
    n_out_files = 1
    
    time.sleep(ARGS.sleep)
    for (i, sub_dir) in enumerate(sub_dirs):
        logging.debug('waiting for %s' % sub_dir)
        with cd(sub_dir):
            # read the end of the file test.out
            out_file = 'par_crysol_%02d.out' % (i+1)
            while tail(out_file, 1) != ['GCRYSOL IS DONE\n']:
                time.sleep(ARGS.sleep)

            print 'sub_run %d is complete' % i
            sub_files = glob.glob('par/crysol/*int')
            error = sub_files.sort()
            for (j, sub_file) in enumerate(sub_files):
                logging.debug('moving %s' % sub_file)
                file_name = sub_file[:-4]
                new_name = ARGS.runname + '_' + str(n_out_files + j).zfill(5)
                logging.debug('moving %s%s to %s/crysol/%s' % (sub_dir, 
                                    file_name, ARGS.runname, new_name))
                os.system('mv %s.int ../crysol/%s.int' % (file_name, new_name))
                os.system('mv %s.log ../crysol/%s.log' % (file_name, new_name))                
                
            os.system('mv *out *dcd *pdb ../')
            n_out_files += len(sub_files)
            
        shutil.rmtree(sub_dir)

def iterate_crysol(sub_dirs, dcd_file_names):
    if os.path.exists('/share/apps/bin/python'):
        #gibbs
        python = '/share/apps/bin/python'
    elif os.path.exists('/usr/bin/python'):
        #my machines
        python = '/usr/bin/python'
    else:
        #others
        python = '/usr/bin/env python'
        
    for (i, sub_dir) in enumerate(sub_dirs):
        run_str = '%s ../../my_crysol_driver.py -r par -pp ./ -p %s -dp ./ -d %s > par_crysol_%02d.out  &' % (python, ARGS.pdb, dcd_file_names[i],i+1)
        # print '\n' + run_str
        with cd(sub_dir):
            os.system(run_str)

def split_dcd():
    import sassie.sasmol.sasmol as sasmol
    
    if os.path.exists(ARGS.runname):
        print 'run folder exists: %s' % ARGS.runname
        print 'select one of the following (0/1/2): quit / move / replace'
        folder = folder_exists()
        result = folder.cmdloop()

    print 'created new run folder: %s' % ARGS.runname
    os.mkdir(ARGS.runname)
    
    mol = sasmol.SasMol(0)
    mol.read_pdb(ARGS.pdb)
    
    mol.read_dcd(ARGS.dcd)
    total_frames = mol.number_of_frames()
    if ARGS.ncpu < 0:
        print 'ncpu: %d < 0,   using |%d| = %d instead' % (ARGS.ncpu, 
                                                ARGS.ncpu, abs(ARGS.ncpu) )
        ARGS.ncpu = abs(ARGS.ncpu)
    n_frames_sub = total_frames/ARGS.ncpu
    last_frame = 0
    sub_dirs = []
    dcd_file_names = []
    for i in xrange(ARGS.ncpu):
        i+=1
        sub_dir = ARGS.runname + '/sub' + str(i).zfill(2) + '/'
        sub_dirs.append(sub_dir)
        os.mkdir(sub_dir)
        os.system('cp %s %s' % (ARGS.pdb, sub_dir))
        with cd(sub_dir):
            if i == ARGS.ncpu:
                n_frames_sub = n_frames_sub + total_frames % ARGS.ncpu
            dcd_file_name = ARGS.dcd[:-4] + '_sub' + str(i).zfill(2) + '.dcd'
            dcd_file_names.append(dcd_file_name)
            first = last_frame
            last = last_frame + n_frames_sub
            mol.write_dcd_frames(dcd_file_name, first, last)
            
        last_frame += n_frames_sub
    print 
    return sub_dirs, dcd_file_names

def main():
    
    #break dcd into N dcds with a folder for each
    sub_dirs, dcd_file_names = split_dcd()
    
    #run crysol instance on each folder
    iterate_crysol(sub_dirs, dcd_file_names)
    
    #collect the results
    collect_crysol(sub_dirs)

    print 'finished all crysol calculations\n     \m/ >.< \m/ '

if __name__ == '__main__':
    import argparse
    
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    ARGS = parse()    
    main()