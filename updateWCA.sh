#!/bin/bash
# $Id: updateWCA.sh,v 1.2 2014-05-23 20:00:27 schowell Exp $
cd /home/schowell/sassie_0.99_rev_1176/sassie/simulate/energy/extensions/non_bonding
cp /home/schowell/Dropbox/gw_phd/code/pylib/sassie/electrostatics.f ./
sudo python setup_electrostatics.py build
sudo chown schowell:schowell -R build/
cp build/lib.linux-i686-2.6/electrostatics.so /home/schowell/Dropbox/gw_phd/code/pylib/sassie/

