#!/bin/bash

# basic script to run the ECE4 job

expname=TEST
ECEDIR=BASEDIR
SRCDIR=$ECEDIR/scripts/platforms
platform=ecmwf-hpc2020-intel+openmpi

se user-config.yml ${expname}.yml $SRCDIR/$platform.yml scriptlib/main.yml --loglevel debug

