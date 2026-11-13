#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to add latitude and longitude to dataset with regular grid 360x180 (x,y)

Author: Alessandro Sozza (CNR-ISAC)
Date: July 2025
"""

# in CDO:
# cdo setgrid,r360x180 woa13-levitus_L31.nc woa13-levitus_L31_r360x180.nc

# example:
# ./add_lonlat.py -f /ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial/woa13-levitus_L31.nc

import sys
import os
import argparse
import xarray as xr
import numpy as np
from scipy.interpolate import interp1d

def main(path):

    # load dataset
    ds = xr.open_dataset(path)

    # create 1d array
    nx, ny = 360, 180
    lon_1d = np.linspace(-179.5, 179.5, nx)
    lat_1d = np.linspace(-89.5, 89.5, ny)

    # compose 2d array
    lon2d, lat2d = np.meshgrid(lon_1d, lat_1d)

    # create variables in dataset
    ds['lon'] = (('y', 'x'), lon2d)
    ds['lat'] = (('y', 'x'), lat2d)

    # add attributes
    ds['lon'].attrs = {'units': 'degrees_east', 'standard_name': 'longitude'}
    ds['lat'].attrs = {'units': 'degrees_north', 'standard_name': 'latitude'}

    # new path
    folder = os.path.dirname(path)
    basename = os.path.basename(path)
    name, ext = os.path.splitext(basename)
    new_path = os.path.join(folder, name + "_r360x180" + ext)

    # write output
    ds.to_netcdf(new_path)
                

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add longitude & latitude to a NetCDF file.")
    parser.add_argument("--file", "-f", required=True, help="NetCDF file")

    args = parser.parse_args()

    main(args.file)

