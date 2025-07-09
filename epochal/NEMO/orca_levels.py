#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to interpolate vertical levels in NEMO datasets.

Author: Alessandro Sozza (CNR-ISAC)
Date: July 2025
"""

# in CDO (not working yet):
# cdo setzaxis,L75.txt woa13-levitus.nc woa13-levitus_L75.nc
# cdo intlevel,zdescription=L31.txt woa13-levitus_L75.nc woa13-levitus_L31.nc
# cdo setgrid,r360x180 woa13-levitus_L31.nc woa13-levitus_L31_r360x180.nc

# example:
# ./orca_levels.py -i /ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial/woa13-levitus.nc -s /ec/res4/hpcperm/itas/data/ece-4-database/nemo/domain/eORCA1/domain_cfg.nc -d /ec/res4/hpcperm/itas/data/ece-4-database/nemo/domain/ORCA2/domain_cfg.nc -o /ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial/woa13-levitus_L31.nc


import sys
import os
import argparse
import xarray as xr
import numpy as np
from scipy.interpolate import interp1d

axis_candidates = {
    'time': ['t', 'time', 'time_counter'],
    'x': ['x', 'lon', 'x_grid_T', 'x_grid_U', 'x_grid_V', 'x_grid_W'],
    'y': ['y', 'lat', 'y_grid_T', 'y_grid_U', 'y_grid_V', 'y_grid_W'],
    'z': ['z', 'lev', 'nav_lev', 'depth', 'deptht', 'depthu', 'depthv', 'depthw']
}

def detect_axis(ds, axis_type, where='dims'):
    
    candidates = axis_candidates.get(axis_type, [])
    if where in ['dims', 'coords']:
        search_space = getattr(ds, where, {})  # ds.dims or ds.coords
    else:        
        search_space = ds.data_vars
        
    for candidate in candidates:
        if candidate in search_space:
            return candidate
    print(f"No {axis_type} {where} found among {candidates}")

    return None

def interpolate(data, old_depths, new_depths, axis):

    new_data = interp1d(old_depths, data, axis=axis, bounds_error=False, fill_value="extrapolate")
    
    return new_data(new_depths)


def main(input_nc, srcdomain_nc, dstdomain_nc, output_nc):
    
    ds = xr.open_dataset(input_nc)
    srcdomain = xr.open_dataset(srcdomain_nc)
    dstdomain = xr.open_dataset(dstdomain_nc)
    
    # assign axis
    axes = {}
    for ax in ['time', 'x', 'y', 'z']:
        axes[ax] = detect_axis(ds, ax, where='dims')
        
    # assign depths
    depth = detect_axis(srcdomain, 'z', where='vars')
    old_depths = srcdomain[depth].values
    new_depths = dstdomain[depth].values
    
    output_vars = {}
    encoding = {}
    for varname in ds.data_vars:
        var = ds[varname]
        
        if axes['z'] not in var.dims:
            continue  # Skip variables without vertical dimension

        if axes['time']:
            data = []
            for t in range(var.sizes[axes['time']]):                
                slice_t = var.isel({axes['time']: t}).values
                interp_slice = interpolate(slice_t, old_depths, new_depths, axis=0)
                data.append(interp_slice)
            new_array = np.stack(data, axis=0)
            new_dims = (axes['time'], axes['z'], axes['y'], axes['x'])
            new_coords = {axes['time']: ds[axes['time']]}
        else:
            data = var.values
            new_array = interpolate(data, old_depths, new_depths, axis=0)
            new_dims = (axes['z'], axes['y'], axes['x'])
            new_coords = {}
            
        output_vars[varname] = xr.DataArray(new_array, dims=new_dims, coords=new_coords, name=varname, attrs=var.attrs)
        
        encoding[varname] = {
            '_FillValue': 9.96921e+36,
            'missing_value': 9.96921e+36,
            'zlib': True,
            'complevel': 4,
            'dtype': 'float32'
        }

    new_ds = xr.Dataset(output_vars)    
    new_ds.attrs = ds.attrs
    
    if axes['time']:
        new_ds.to_netcdf(output_nc, encoding=encoding, unlimited_dims=axes['time'])
    else:
        new_ds.to_netcdf(output_nc, encoding=encoding)        
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Interpolate vertical levels using two domain_cfg.nc files.")
    parser.add_argument("--infile", "-i", required=True, help="Input NetCDF file to interpolate")
    parser.add_argument("--srcdomain", "-s", required=True, help="Source domain_cfg.nc with original vertical levels")
    parser.add_argument("--dstdomain", "-d", required=True, help="Destination domain_cfg.nc with target vertical levels")
    parser.add_argument("--outfile", "-o", required=True, help="Output NetCDF file")

    args = parser.parse_args()

    main(args.infile, args.srcdomain, args.dstdomain, args.outfile)

