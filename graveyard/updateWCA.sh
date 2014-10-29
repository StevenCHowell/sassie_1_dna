#!/bin/bash
#
# Author:  Steven C. Howell
# Purpose: build the WCA energy calculator
# Created: January 2014
#
# $Id$
#
cd /home/schowell/sassie_0.99_rev_1176/sassie/simulate/energy/extensions/non_bonding
cp /home/schowell/Dropbox/gw_phd/code/pylib/sassie/electrostatics.f ./
sudo python setup_electrostatics.py build
sudo chown schowell:schowell -R build/
cp build/lib.linux-i686-2.6/electrostatics.so /home/schowell/Dropbox/gw_phd/code/pylib/sassie/

