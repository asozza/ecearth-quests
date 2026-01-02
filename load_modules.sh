#!/bin/bash

module load python3/3.12
module load cdo nco ncview
module load intel/2021.4.0 intel-mkl/19.0.5 prgenv/intel hdf5 netcdf4

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/apps/netcdf4/4.9.1/INTEL/2021.4/lib:/usr/local/apps/hdf5/1.12.2/INTEL/2021.4/lib

