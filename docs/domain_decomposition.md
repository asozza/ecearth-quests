
## Best domain decomposition

The aim is to print a list of possible MPI tasks, ensuring best domain decompositions. 
We need `nemo.exe` to run NEMO standalone and a domain configuration file `domain_cfg.nc`.

Set in `&nammpp` `ln_listonly = .true.`:

```
!-----------------------------------------------------------------------
&nammpp        !   Massively Parallel Processing
!-----------------------------------------------------------------------
   ln_listonly =  .true.   !  do nothing else than listing the best domain decompositions (with land domains suppression)
   !                       !  if T: the largest number of cores tested is defined by max(mppsize, jpni*jpnj)
   ln_nnogather =  .true.  !  activate code to avoid mpi_allgather use at the northfold
   jpni        =   0       !  number of processors following i (set automatically if < 1), see also ln_listonly = T
   jpnj        =   0       !  number of processors following j (set automatically if < 1), see also ln_listonly = T
   nn_hls      =   1       !  halo width (applies to both rows and columns)
   nn_comm     =   1       !  comm choice
/
```

Here's the sbatch script to run NEMO in parallel:

```
#!/bin/bash

#SBATCH --job-name=nemo
#SBATCH --output=log_nemo_%j.out
#SBATCH --error=log_nemo_%j.err    
#SBATCH --qos=np
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=01:00:00
#SBATCH --account=hpc-account

# set environment
manager=micromamba
env=ece4 
file=mamba

basepath=$(manager env list | grep base | awk '{print $2}')
source $basepath/etc/profile.d/$file.sh

# load environment
micromamba activate ece4

# load modules
module reset
module load prgenv/intel
module load intel/2021.4.0
module load hpcx-openmpi/2.9.0
module load intel-mkl/19.0.5
module load fftw/3.3.9
module load netcdf4-parallel/4.9.1
module load cmake/3.20.2

# add libraries to path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/apps/netcdf4-parallel/4.9.1/INTEL/2021.4/HPCX/2.9/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ec/res4/hpcperm/itas/src/gitlab/ecearth4-fork/sources/oasis3-mct-5.2/arch_ecearth/lib

# in the current directory
cd .

# run nemo stand-alone
mpirun ./nemo.exe
```
