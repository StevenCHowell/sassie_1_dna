#!/bin/bash
# $Id: updateCollision.sh,v 1.3 2014-08-22 18:05:12 schowell Exp $
ORIGINAL=`pwd`
# DIR="../simulate/energy/extensions/non_bonding/"
DIR="/usr/local/lib/python2.7/dist-packages/sassie/simulate/energy/extensions/non_bonding/"
cp collision.f $DIR -f
cp setup_collision.py $DIR -f
cd $DIR
sudo python setup_collision.py build
# sudo chown schowell:schowell -R build/
sudo chown root:root -R build/
cp build/lib.linux-x86_64-2.7/collision.so $ORIGINAL
cd $ORIGINAL
