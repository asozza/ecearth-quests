#!/usr/bin/env python
"""
Function to create corner coordinates for a grid based on input ICMGG files from OIFS.
Based on the original script from Florian Ziemen.

Paolo Davini, Apr 2024
"""

import os
import netCDF4 as nc
import numpy as np
from utils import extract_grid_info, spectral2gaussian
import cdo
cdo = cdo.Cdo()
#cdo.debug = True


resolutions = ["TL63L31", "TL159L91"]
oifs_dir = "/lus/h2resw01/hpcperm/ccpd/ECE4-DATA/oifs"
tgt_dir = "/ec/res4/scratch/itmn/IFS-masked"


for resolution in resolutions:
    print('Processing resolution:', resolution)

    kind, spectral, vertical =  extract_grid_info(resolution)
    infile_name = f"{oifs_dir}/{resolution}/19900101/ICMGGECE4INIT"
    netcdf_name = f"{tgt_dir}/{resolution}-tmp.nc"
    gaussian_name = f"{tgt_dir}/{resolution}-gaussian.nc"
    outfile_name = f"{tgt_dir}/T{kind}{spectral}_grid.nc"
    outfile_masked_name = f"{tgt_dir}/T{kind}{spectral}_grid_masked.nc"
    variable_name = "var172" #land-sea mask, but anything else will work
    gaussian = spectral2gaussian(spectral, kind)

    # convert to netcdf
    print(f"Converting GRIB {infile_name} to NetCDF {netcdf_name}")
    cdo.selname(variable_name, input=infile_name, output=netcdf_name, options="-f nc4")

    # there is a strange bug in cdo. latitude are not recognized in the original grid file
    # however, we know them from gaussian grid associated
    print(f"Creating gaussian grid {gaussian_name}")
    cdo.const(f"1,n{gaussian}", output=gaussian_name, options="-f nc4")
    #ubprocess.call(["cdo", "-f", "nc4", f"const,1,n{gaussian}", gaussian_name])

    # load the latitudes from the new file
    infile = nc.Dataset(gaussian_name)
    variables = infile.variables
    lat = variables["lat"]

    # load netcdf
    infile = nc.Dataset(netcdf_name)
    variables = infile.variables

    print(variables[variable_name].shape)
    total_points=variables[variable_name].shape[-1]
    rp = variables["reduced_points"][:]

    if lat.shape[0] != len(rp):
        raise ValueError("Number of latitudes does not match number of reduced points")


    lats=np.zeros(total_points)
    lons=np.zeros(total_points)
    lons_left=np.zeros(total_points)
    lons_right=np.zeros(total_points)
    lats_upper=np.zeros(total_points)
    lats_lower=np.zeros(total_points)

    # the hypothesis is that the latitudan bands are equally spaced, so that corners lies on the midpoints
    print("Creating corner coordinates...")
    lat_upper=lat[:]
    lat_upper[:-1]=lat_upper[:-1]+.5*(lat_upper[:-1]-lat_upper[1:])
    lat_upper[-1]=-lat_upper[1]
    lat_lower=lat[:]
    lat_lower[:-1]=lat_upper[1:]
    lat_lower[-1]=-lat_upper[0]

    covered=0
    # similar assumpution is done for the longitudes
    for j, num_points in enumerate(rp):
        lats[covered:covered+num_points] = lat[j]
        lons[covered:covered+num_points] = np.arange(num_points)/num_points*360
        lons_left[covered:covered+num_points] = (np.arange(num_points)-.5)/num_points*360
        lons_right[covered:covered+num_points] = (np.arange(num_points)+.5)/num_points*360
        lats_upper[covered:covered+num_points] = lat_upper[j]
        lats_lower[covered:covered+num_points] = lat_lower[j]
        covered += num_points

    for x in [ lons, lons_left, lons_right]:
        x[x>180] = x[x>180]-360.


    corners_lat=np.array([lats_upper,lats_upper,lats_lower,lats_lower]).transpose()
    corners_lon=np.array([lons_left,lons_right, lons_right, lons_left]).transpose()

    print("Writing output file...", outfile_name)
    ds = nc.Dataset(outfile_name, "w", format="NETCDF4")
    cells = ds.createDimension("rgrid",total_points)
    nv = ds.createDimension("nv", 4)
    clon = ds.createVariable("clon","f8", ("rgrid",))
    clat = ds.createVariable("clat","f8", ("rgrid",))
    lml = ds.createVariable("lml","f4", ("rgrid",)) # fake variable for CDI
    clon_bnds = ds.createVariable("clon_bnds","f8", ("rgrid", "nv"))
    clat_bnds = ds.createVariable("clat_bnds","f8", ("rgrid", "nv"))
    for x in [ clon, clat, lml]:
        x.units="radian"
        x.coordinates="clat clon"
    clon.bounds = "clon_bnds"
    clat.bounds = "clat_bnds"
    clon.standard_name = "longitude"
    clat.standard_name = "latitude"

    clon[:]=lons*np.pi/180.
    clat[:]=lats*np.pi/180.
    lml[:]=clat[:]*clon[:]
    clon_bnds[:]=corners_lon*np.pi/180
    clat_bnds[:]=corners_lat*np.pi/180
    ds.close()

    print("Writing masked output file...", outfile_masked_name)
    ds_masked = nc.Dataset(outfile_masked_name, "w", format="NETCDF4")
    cells_masked = ds_masked.createDimension("rgrid",total_points)
    nv_masked = ds_masked.createDimension("nv", 4)
    clon_masked = ds_masked.createVariable("clon","f8", ("rgrid",))
    clat_masked = ds_masked.createVariable("clat","f8", ("rgrid",))
    variable_mask = variables[variable_name][:]
    mask = ds_masked.createVariable("lsm","f4", ("rgrid",))
    mask[:] = variable_mask[:]
    clon_bnds_masked = ds_masked.createVariable("clon_bnds","f8", ("rgrid", "nv"))
    clat_bnds_masked = ds_masked.createVariable("clat_bnds","f8", ("rgrid", "nv"))
    for x in [ clon_masked, clat_masked, mask]:
        x.units="radian"
        x.coordinates="clat clon"
    clon_masked.bounds = "clon_bnds"
    clat_masked.bounds = "clat_bnds"
    clon_masked.standard_name = "longitude"
    clat_masked.standard_name = "latitude"

    clon_masked[:]=lons*np.pi/180.
    clat_masked[:]=lats*np.pi/180.
    mask[:]=variables[variable_name][:]
    clon_bnds_masked[:]=corners_lon*np.pi/180
    clat_bnds_masked[:]=corners_lat*np.pi/180
    ds_masked.close()

    print("Cleaning up...")
    os.remove(netcdf_name)
    os.remove(gaussian_name)
