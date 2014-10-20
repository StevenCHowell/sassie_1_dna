#!/bin/bash
#
# Author:  Steven C. Howell
# Purpose: build a local version of the coulomb calculator
# Created: August 2014
#
# $Id: updateCollisionLocal.sh 77 2014-10-17 16:19:32Z schowell $

/usr/bin/env python setup_coulomb.py build
cp build/*/coulomb.so ./

if [ ! -f ./coulomb.so ] ; then
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! failed to build 'coulomb.so' using /usr/bin/env python !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! retrying using /share/apps/bin/python                    !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
    sleep 1

    /share/apps/bin/python setup_coulomb.py build
    cp build/*/coulomb.so ./

    if [ ! -f ./coulomb.so ] ; then
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! still could not identify where the executable is located !!!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! compilation likely failed both times                     !!!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! Please inspect compilation output for bugs               !!!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! if 'coulomb.so' exists, please move it to ./           !!!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
    else
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! sucessfully created 'coulomb.so' and moved it to ./ !!!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
    fi
else
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! sucessfully created 'coulomb.so' and moved it to ./ !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
fi

#s USER=`whoami`
#s if [ "$USER" == schowell ] ; then
#s     rm -rf build
#s fi
