#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: Run ds_dna_monte_carlo.py using gui-like input
# Created: 30 Octomber 2014
#
# $Id$
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import dna.ds_dna_monte_carlo as ddmc
import dna.special.input_filter as input_filter

def prepare_dna_mc_input(variables):
    # ~~~~Generate 'flex_resid' list~~~~
    theta_max             = variables['theta_max'][0]
    theta_z_max           = variables['theta_z_max'][0]
    n_flex_regions        = variables['n_flex_regions'][0]
    first_res_per_region  = variables['first_res_per_region'][0]
    n_cont_res_per_region = variables['n_cont_res_per_region'][0]

    assert len(theta_max) == n_flex_regions, 'theta_max should have %d values' % n_flex_regions
    assert len(theta_z_max) == n_flex_regions, 'theta_z_max should have %d values' % n_flex_regions
    assert len(first_res_per_region) == n_flex_regions, 'first_res_per_region should have %d values' % n_flex_regions
    assert len(n_cont_res_per_region) == n_flex_regions, 'n_cont_res_per_region should have %d values' % n_flex_regions
    
    flex_resids = []
    for i in xrange(n_flex_regions):
        first = first_res_per_region[i]
        last  = first_res_per_region[i] + n_cont_res_per_region[i]
        flex_resids.append(range(first, last))

    variables['flex_resids'] = (flex_resids, 'list_of_lists')

    # ~~~~Generate 'pro_groups' list~~~~
    n_pro_groups = variables['n_pro_groups'][0]
    assert n_pro_groups <= n_flex_regions, '# of protein groups must be <= to # of flexible DNA regions'
    pro_groups = []
    for i in xrange(n_pro_groups):
        pro_groups.append(variables['pro_group'+str(i+1)][0])
    variables['pro_groups'] = (pro_groups, 'list_of_lists')
        
    # ~~~~Generate 'dna_resids' list~~~~
    dna_resids = [variables['dna1_resids'][0], variables['dna2_resids'][0]]
    variables['dna_resids'] = (dna_resids, 'list_of_lists')
    
    return variables

if __name__ == "__main__":

    svariables = {}

    # User Input
    svariables['outfile']   = ('new_dsDNA.dcd', 'string')
    svariables['pdbfile']   = ('new_dsDNA.pdb', 'string')
    svariables['trials']    = ('100',           'int')
    svariables['goback']    = ('50',            'int')
    
    # Molecule Specific Input
    svariables['n_flex_regions']        = ('3', 'int')
    svariables['theta_max']             = ('25, 10, 5', 'float_array')
    svariables['first_res_per_region']  = ('2, 10, 20', 'int_array')
    svariables['n_cont_res_per_region'] = ('2, 5, 7', 'int_array')
    svariables['align_low_res']         = ('16', 'int')  # not yet implemented
    svariables['align_high_res']        = ('19', 'int')  # not yet implemented
    # ~~~ TO BE ADDED TO THE GUI ~~~ # 
    svariables['dna_segnames']          = ('DNA1, DNA2', 'string_array')  
    svariables['dna1_resids']           = ('1, 30', 'int_array')
    svariables['dna2_resids']           = ('30, 1', 'int_array')
    
    # this is how protein input should be for an example case of 2 protein groups, 3 flexible DNA groups
    svariables['n_pro_groups']          = ('3', 'int')
    # 1- after the first flexible DNA region
    svariables['pro_group1']            = ('A0, B0, C0, D0, E0, F0, G0, H0', 'string_array')
    # place holder
    svariables['pro_group2']            = ('', 'string_array')
    # 2 -after the second flexible DNA region
    svariables['pro_group3']            = ('A1, B1, C1, D1, E1, F1, G1, H1', 'string_array')

    # when there is a protein group before any flexible regions, 
    # a flex region with 0 continuous residues must be entered (the first residue can be arbitrary), i.e.:
    # > svariables['first_res_per_region']  = ('1, 10, 20', 'int_array')
    # > svariables['n_cont_res_per_region'] = ('0, 5, 7', 'int_array')
    
    # Non-Standard/Specialized Inputs (normal users should use a default value)
    if 'theta_z_max' in locals():
        s_theta_z_max = ''
        for (i, theta) in enumerate(theta_z_max):
            if i > 0:
                s_theta_z_max += ', '
            s_theta_z_max += str(theta)    
    else:
        s_theta_z_max = svariables['theta_max'][0]

    svariables['theta_z_max']   = (s_theta_z_max, 'float_array') # provide 's_theta_z_max' to use a different theta max for twisting vs bending

    svariables['bp_per_bead']   = ('1',     'int') # set to N > 1 to have N base-pairs coarse-grained into 1 bead
    svariables['softrotation']  = ('1',     'int') # set to N > 1 to apply rotations averaged over N coars-grained beads
    svariables['n_dcd_write']   = ('1',     'int') # set to N > 1 to only record structures after every N trials
    svariables['seed']          = ('0',     'int') # goback random seed
    svariables['temperature']   = ('300',   'float') # may be incorperated for electrostatics
    svariables['debug']         = ('False', 'bool') # 'True' will display debug output
    svariables['write_flex']    = ('False', 'bool') # 'True' will generate a file labeling paird DNA base pairs for all flexible beads (useful for calculating dihedral angles)
    svariables['keep_cg_files'] = ('False', 'bool') # 'True' will keep the coarse-grain DNA and protein pdb and dcd files
    svariables['keep_unique']   = ('True',  'bool') # 'False' will keep duplicate structure when move fails
    svariables['rm_pkl']        = ('False', 'bool') # 'True' will remove the coarse-grain pkl file forcing a re-coarse-graining (can be time consuming often necessary)

    # Parse input
    error, variables = input_filter.type_check_and_convert(svariables)
    variables = prepare_dna_mc_input(variables)

    # Run 
    if 'finished' == ddmc.main(variables):
        trials = variables['trials'][0]
        print '\nFinished %d successful DNA moves! \n\m/ >.< \m/' % trials
    else:
        print '\n~~~ Did not run sucessfully. ~~~'
        print '----> Check input parameters. <----'
