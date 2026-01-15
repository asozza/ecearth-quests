#!/bin/bash

# tentative full workflow to create PALEORCA input data from raw data files

# reset modules and load required ones
module reset
module load prgenv/intel intel/2021.4.0 intel-mkl/19.0.5 hpcx-openmpi/2.9.0
module load hdf5-parallel/1.12.2 netcdf4-parallel/4.9.1 ecmwf-toolbox/2023.04.1.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF4_PARALLEL_DIR/lib:$ECCODES_DIR/lib:$HDF5_DIR/lib:$HPCPERM/ecearth4/revisions/main/sources/oasis3-mct-5.2/arch_ecearth/lib

OUTPUTDIR=/lus/h2resw01/hpcperm/ccpd/EPOCHAL_v2
mkdir -p $OUTPUTDIR

# Create halo for coordinates
#cdo sethalo,-1,-1 coords_ori.nc coords_halo.nc

# run the script to create the bounds: grid staggering is controversial, using T grid here
COORDSDIR=/lus/h2resw01/hpcperm/ccpd/EPOCHAL/PALEORCA-coords
staggering=T
python3 orca_bounds_new.py --stagg $staggering --no-level $COORDSDIR/coords_halo.nc $COORDSDIR/coords_bounds_$staggering.nc

# remap bathymetry from eORCA1 to PALEORCA grid using nearest neighbor
source_grid=eORCA1
DOMAINDIR=/lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain
GRIDDIR=/lus/h2resw01/hpcperm/ccpd/EPOCHAL/ORCA
cdo remapnn,${COORDSDIR}/coords_bounds_$staggering.nc -setgrid,$GRIDDIR/${source_grid}/${source_grid}_grid_${staggering}.nc -selname,bathy_metry,nav_lon,nav_lat $DOMAINDIR/${source_grid}/domain_cfg.nc $OUTPUTDIR/PALEORCA_bathy_metry_from_${source_grid}.nc

# remape the Herold et al. (2014) Eocene topography to PALEORCA grid using nearest neighbor
HEROLDDIR=/lus/h2resw01/hpcperm/ccpd/EPOCHAL/Herold_etal_2014
cdo chname,topo,bathy_metry -mulc,-1 -setrtoc,0,10000,0 -remapnn,${COORDSDIR}/coords_bounds_$staggering.nc -selname,topo,lon,lat ${HEROLDDIR}/herold_etal_eocene_topo_1x1.nc $OUTPUTDIR/PALEORCA_bathy_metry_fromHerold.nc

