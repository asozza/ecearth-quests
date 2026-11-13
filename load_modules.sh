#!/bin/bash

# Working ONLY on ATOS HPC2020

module load prgenv/intel
module load intel/2021.4.0
module load hpcx-openmpi/2.9.0
module load intel-mkl/19.0.5
module load fftw/3.3.9
module load cmake/3.20.2
module load netcdf4-parallel/4.9.1
module load hdf5-parallel/1.10.6
module load python3/3.12
module load cdo

# needed for running REBUILD_NEMO
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF4_PARALLEL_DIR/lib:$ECCODES_DIR/lib:$HDF5_DIR/lib
