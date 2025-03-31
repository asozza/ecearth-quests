#!/bin/bash

# before creating a new folder check paths in user-config.yml
# remember also that ini_dir is defined in /sources/se/ecmwf

jobname=$1
kind=$2

# only this one is supported
machine=ecmwf-hpc2020-intel+openmpi

# please define where the source code is
#ecedir=$HPCPERM/src/ecearth4-epochal
ecedir=$HPCPERM/src/gitlab/ecearth4-fork
#ecedir=$HPCPERM/ecearth4/revisions/main

# please define where the jobs are
expdir=$HPCPERM/ecearth4/jobs
default=$ecedir/scripts/runtime
rundir=$SCRATCH/ece4

# hard-coded to use Ale new updated files
inidir=/ec/res4/hpcperm/itas/data/ECE4-DATA


if [ -z $jobname ] ; then
	echo " usage ./create_new_job.sh jobname kind ('amip' or 'cpld') "
	exit 1
fi

if [ -z $kind ] ; then
        echo " usage ./create_new_job.sh jobname kind ('amip' or 'cpld') "
        exit 1
fi


mkdir -p $expdir/$jobname
cp -r $default/scriptlib $expdir/$jobname
cp -r $default/templates $expdir/$jobname

cp experiment-config-$kind.yml $expdir/$jobname/$jobname.yml
sed -i "s/TEST/${jobname}/g" $expdir/$jobname/$jobname.yml

cp user-config-example.yml $expdir/$jobname/user-config.yml
sed -i "s@RUNDIR@${rundir}@g" $expdir/$jobname/user-config.yml
sed -i "s@BASEDIR@${ecedir}@g" $expdir/$jobname/user-config.yml
sed -i "s@INIDIR@${inidir}@g" $expdir/$jobname/user-config.yml

cp launch.sh $expdir/$jobname
sed -i "s@TEST@${jobname}@g" $expdir/$jobname/launch.sh
sed -i "s@BASEDIR@${ecedir}@g" $expdir/$jobname/launch.sh

