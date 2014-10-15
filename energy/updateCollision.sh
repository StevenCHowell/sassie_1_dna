#!/bin/bash
#
# Author:  Steven C. Howell
# Purpose: build the calculator for calculating energy between each pair of atoms
# Created: May 2014
# $Id$
#

ORIGINAL=`pwd`
DIR="/usr/local/lib/python2.7/dist-packages/sassie/simulate/energy/extensions/non_bonding/"
cp collision.f $DIR -f
cp setup_collision.py $DIR -f
cd $DIR
sudo python setup_collision.py build
sudo chown root:root -R build/
cp build/*/collision.so $ORIGINAL
cd $ORIGINAL
