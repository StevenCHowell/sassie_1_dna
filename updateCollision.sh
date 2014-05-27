#!/bin/bash
# $Id: updateCollision.sh,v 1.2 2014-05-27 18:03:08 schowell Exp $
ORIGINAL=`pwd`
DIR="../simulate/energy/extensions/non_bonding/"
cp collision.f $DIR -f
cp setup_collision.py $DIR -f
cd $DIR
sudo python setup_collision.py build
sudo chown schowell:schowell -R build/
cp build/lib.linux-x86_64-2.7/collision.so $ORIGINAL
cd $ORIGINAL
