#!/bin/bash
#
# Author:  Steven C. Howell
# Purpose: build a local version of the collision calculator
# Created: August 2014
#
# $Id$

/usr/bin/env python setup_collision.py build

cp build/*/collision.so ../

if [ ! -f ../collision.so ] ; then
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! failed to identify where the executable is located !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! compilation failed or your system is unfamiliar    !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!!! Please find and move 'collision.so' to ../ folder  !!!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
fi

#s USER=`whoami`
#s if [ "$USER" == schowell ] ; then
#s     rm -rf build
#s fi
