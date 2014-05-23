#!/bin/bash
# $Id: updateCollision.sh,v 1.1 2014-05-23 20:00:27 schowell Exp $
ORIGINAL=`pwd`
DIR="../simulate/energy/extensions/non_bonding/"
cp collision.f $DIR -f
cp setup_collision.py $DIR -f
cd $DIR
sudo python setup_collision.py build
sudo chown schowell:schowell -R build/
cp build/lib.linux-i686-2.6/electrostatics.so $ORIGINAL

