#!/bin/bash
# $Id: updateCollision.sh,v 1.4 2014-10-14 15:45:38 schowell Exp $
ORIGINAL=`pwd`
DIR="/usr/local/lib/python2.7/dist-packages/sassie/simulate/energy/extensions/non_bonding/"
cp collision.f $DIR -f
cp setup_collision.py $DIR -f
cd $DIR
sudo python setup_collision.py build
sudo chown root:root -R build/
cp build/*/collision.so $ORIGINAL
cd $ORIGINAL
