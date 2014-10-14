#!/bin/bash
#
# Author:  Steven C. Howell
# Purpose: update the files in trunk version of SASSIE to allow for calculation dihedral angles of dsDNA
# Created: September 2014
#
# $Id$
#
cp mask.c /home/schowell/data/myPrograms/svn_utk/trunk/sassie/sasmol/extensions/mask/mask.c
cp dihedral_monte_carlo.py /home/schowell/data/myPrograms/svn_utk/trunk/sassie/simulate/monte_carlo/monomer/dihedral_monte_carlo.py
cp dihedral_rotate.py /home/schowell/data/myPrograms/svn_utk/trunk/sassie/simulate/monte_carlo/monomer/dihedral_rotate.py
