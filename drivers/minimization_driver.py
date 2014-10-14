#!/usr/bin/env python
# Author:  Steven C. Howell
# Purpose: script for running SASSIE minimizer with command line input
# Created: 9 September 2014                                                   
# $Id$

import os, sys, shutil, glob
import multiprocessing
import sassie.interface.input_filter as input_filter

import sassie.interface.minimize_filter as minimize_filter
import sassie.simulate.namd.namd_minimize as minimize

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'Minimize Structures'
        #epilog = 'no epilog found'
    )
    # standard input
    parser.add_argument("-r", "--runname", type=str, help = ("folder to store the minimization results"))
    parser.add_argument("-i", "--infile", type=str, help =("dcd/pdb file that contains the structure/s to minimize"))
    parser.add_argument("-p", "--pdbfile", type=str, help = ("pdb file containing the structure info (enables using dcd files)"))
    parser.add_argument("-o", "--outfile", type=str, help = ("name for output dcd file containing all the minimized structures"))
    parser.add_argument("-psf", "--psffile", type=str, help = ("charmm psf file containing force info for molecular bonds"))

    # default usually sufficient
    parser.add_argument("--parmfile", type=str, default = '/usr/local/bin/toppar/par_all27_prot_na.inp', help = ("topology file containing molecule info"))
    parser.add_argument("--nsteps", type=str, default = '2000', help = ("number of minimization steps"))
    parser.add_argument("--ncpu", type=str, default = '8', help = ("number of cpus to use for minimization"))
    parser.add_argument("--keepout", type=str, default = '1', help = ("keep run output"))
    parser.add_argument("--dcdfreq", type=str, default = '20', help = ("frequency that frames are written during the individual minimization (should be an integer devisor of the total # of minimization steps)"))    

    # non-standard / specialized input
    parser.add_argument("--md", type=str, default = '2', help = ("0: min, 1: min-md, 2: min-md-min"))
    parser.add_argument("--mdsteps", type=str, default = '100', help = ("# of md steps"))
    parser.add_argument("--dielect", type=str, default = '80.0', help = ("solvent dielectric"))    
    parser.add_argument("--temperature", type=str, default = '300.0', help = ("temperature"))      
    
    return parser.parse_args()

class Drv():
    module = 'minimization'
    
    def run_me(self):
        import os

        #### BEGIN USER EDIT
        #### BEGIN USER EDIT
        #### BEGIN USER EDIT

        runname = 'dsDNA60bps'
        infile = './140904_094120_new_dsDNA60.dcd'
        pdbfile = './new_dsDNA60.pdb' # this is the atom info file (for use with a dcd)
        outfile = 'dsDNA60_100s_min.dcd'
        psffile = './new_dsDNA60.psf'

        nsteps = '2000'
        parmfile = '/usr/local/bin/toppar/par_all27_prot_na.inp'
        ncpu = '4'
        keepout = '1'
        dcdfreq = '20'
        infiletype = infile[-3:] #'dcd'  # need to change this to dcd if I input a dcd file 

        md = '0'
        mdsteps = '100'
        dielect = '80.0'
        temperature = '300.0'

        #### END USER EDIT
        #### END USER EDIT
        #### END USER EDIT

        uname = os.popen("whoami").read()
        if 'schowell\n' == uname:
            ARGS = parse()           
            runname     = ARGS.runname    
            infile      = ARGS.infile     
            infiletype  = ARGS.infile[-3:]
            pdbfile     = ARGS.pdbfile    
            outfile     = ARGS.outfile    
            nsteps      = ARGS.nsteps     
            parmfile    = ARGS.parmfile   
            psffile     = ARGS.psffile    
            ncpu        = ARGS.ncpu       
            keepout     = ARGS.keepout    
            dcdfreq     = ARGS.dcdfreq    
            md          = ARGS.md         
            mdsteps     = ARGS.mdsteps    
            dielect     = ARGS.dielect    
            temperature = ARGS.temperature

        svariables={}

        svariables['runname']           = (runname,'string')
        svariables['infile']            = (infile,'string')
        svariables['pdbfile']           = (pdbfile,'string')
        svariables['outfile']           = (outfile,'string')
        svariables['nsteps']            = (nsteps,'int')
        svariables['parmfile']          = (parmfile,'string')
        svariables['psffile']           = (psffile,'string')
        svariables['ncpu']              = (ncpu,'int')
        #svariables['energyfile']       = (energyfile,'string')
        svariables['keepout']           = (keepout,'int')
        svariables['dcdfreq']           = (dcdfreq,'int')
        svariables['infiletype']        = ('dcd','string')

        svariables['md']                = (md,'int')
        svariables['mdsteps']           = (mdsteps,'int')
        svariables['dielect']           = (dielect,'float')
        svariables['temperature']       = (temperature,'float')


        error,variables=input_filter.type_check_and_convert(svariables)

        if(len(error) != 0):
            print 'error = ',error
            sys.exit()
        else:
            error=minimize_filter.check_minimize(variables)

            if(len(error) != 0):
                print 'error = ',error
                sys.exit()


        runname=variables['runname'][0]

        import multiprocessing
        import shutil, os
        if os.path.exists(runname+'/'+self.module):
            shutil.rmtree(runname+'/'+self.module)

        txtQueue=multiprocessing.JoinableQueue()
        minimize.minimize(variables,txtQueue)

 

if __name__=='__main__':
    import argparse
    o=Drv()
    o.run_me()
