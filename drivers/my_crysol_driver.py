import sys
import sassie.interface.input_filter as input_filter
import sassie.calculate.gcrysol as gcrysol

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'Minimize Structures'
        #epilog = 'no epilog found'
    )
    # standard input
    parser.add_argument( "-r", "--runname", type=str, help = ("folder to store the minimization results"))
    parser.add_argument( "-d", "--dcdfile", type=str, help = ("dcd/pdb file containing the structures for scattering calculation"))
    parser.add_argument("-dp", "--dcdpath", type=str, default = 'minimization/', help =("path to the dcd/pdb file containing the structures for scattering calculation"))
    parser.add_argument( "-p", "--pdbfile", type=str, default = 'new_c11h5.pdb', help = ("pdb file containing the structure info (enables using dcd files)"))
    parser.add_argument("-pp", "--pdbpath", type=str, default = 'minimization/', help = ("path pdb file containing the structure info (enables using dcd files)"))

    return parser.parse_args()



class Drv():

    module = 'crysol'

    def run_me(self):
        import os

        #### BEGIN USER EDIT
        #### BEGIN USER EDIT
        #### BEGIN USER EDIT

        # unique to run
        runname='test_crysol'
        dcdpath='./'
        dcdfile='test_crysol.dcd'
        pdbpath='./'
        pdbfile='new_dsDNA60.pdb'

        # system defaults
        #cryexe='/share/apps/bin/crysol.exe'
        cryexe='/usr/local/bin/crysol.exe'
        delafs='1'

        # my defaults
        # maxh='5'
        # fib='5'
        maxh='15'
        fib='18'
        maxs='0.199'
        numpoints='200'

        # program defaults
        option='0'
        contrast='0.03'
        edensolv='0.334'
        hydrogens='N'

        #### END USER EDIT
        #### END USER EDIT
        #### END USER EDIT

        uname = os.popen("whoami").read()
        if 'schowell\n' == uname:
            ARGS = parse()
            if ARGS.dcdfile != None and ARGS.runname != None:
                print 'loading parameters from command line'
                runname = ARGS.runname
                dcdpath = ARGS.dcdpath
                dcdfile = ARGS.dcdfile
                pdbpath = ARGS.pdbpath
                pdbfile = ARGS.pdbfile
            else:
                print 'using parameters from driver script'

        svariables={}

        svariables['runname']   = (runname,'string')
        svariables['dcdpath']   = (dcdpath,'string')
        svariables['dcdfile']   = (dcdfile,'string')
        svariables['pdbpath']   = (pdbpath,'string')
        svariables['pdbfile']   = (pdbfile,'string')
        svariables['cryexe']    = (cryexe,'string')
        svariables['delafs']    = (delafs,'int')
        svariables['option']    = (option,'string')
        svariables['maxh']      = (maxh,'string')
        svariables['fib']       = (fib,'string')
        svariables['maxs']      = (maxs,'string')
        svariables['numpoints'] = (numpoints,'string')
        svariables['contrast']  = (contrast,'string')
        svariables['edensolv']  = (edensolv,'string')
        svariables['hydrogens'] = (hydrogens,'string')

        error,self.variables=input_filter.type_check_and_convert(svariables)
        import shutil, os
        if os.path.exists(runname+'/'+self.module):
            shutil.rmtree(runname+'/'+self.module)

        if(len(error)>0):
            print 'error = ',error
            sys.exit()

        runname=self.variables['runname'][0]

        import multiprocessing
        import sassie.tools.center as center

        txtQueue=multiprocessing.JoinableQueue()
        gcrysol.gcrysol(self.variables,txtQueue)

    def verify_me(self):
        module = self.module

        import os, filecmp,glob
        result_path = './'+self.variables['runname'][0]+'/'+module+'/'
        expected_path = './expected_results/'+module+'/'
        result_files = glob.glob(result_path+'*')
        result_files = [os.path.split(result_files[i])[1] for i in range(len(result_files))]
        expected_files = glob.glob(expected_path+'*')
        expected_files = [os.path.split(expected_files[i])[1] for i in range(len(expected_files))]
        print '\nvariables: ',self.variables,'\n\nresult_files:   ',result_path, '\n', result_files, '\n\nexepected_files:',expected_path, '\n',expected_files

        from util import FileCmp
        flag = True
        for ifile in expected_files:
            if ifile in result_files:
                print '\ncomparing ',expected_path+ifile, result_path+ifile,
                if ifile[-4:]=='.log':
                    flag = (flag and FileCmp.cmp_skip(expected_path+ifile, result_path+ifile, [1, 4, 23, 58]))
                else:
                    flag = (flag and filecmp.cmp(expected_path+ifile, result_path+ifile))
                print '\n...to be...',flag
                if flag==False:
                    return False
            else:
                return False
        return flag


if __name__=='__main__':
    import argparse
    o=Drv()
    o.run_me()
    # o.verify_me()
