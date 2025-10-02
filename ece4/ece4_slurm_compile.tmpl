#!/bin/bash

#SBATCH --job-name=compile_ece4
#SBATCH --qos=nf
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=65G
#SBATCH --output=log_compile_ece4.%j.out
#SBATCH --error=log_compile_ece4.%j.out
#SBATCH --account=spitsozz

ECEDIR=/lus/h2resw01/hpcperm/ecme3497/ec-earth-4-fork
platform=ecmwf-hpc2020-intel+openmpi.yml
ECESOURCES=$ECEDIR/scripts/build

# find mamba/conda (to be refined)
whereconda=$(which conda | rev | cut -f 3-10 -d"/" | rev)
source $whereconda/etc/profile.d/conda.sh

conda activate ece4
cd $ECESOURCES

components=""

if [ -z "$components" ]; then
    se user-settings.yml ../platforms/$platform compile-components.yml
else
    for component in $components; do
        se user-settings.yml ../platforms/$platform compile-$component.yml
    done
fi


