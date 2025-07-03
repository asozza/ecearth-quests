#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to interpolate vertical levels in NEMO datasets.

Author: Alessandro Sozza (CNR-ISAC)
Date: July 2025
"""

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
    search_space = getattr(ds, where, {})  # ds.dims or ds.coords

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

    depth = detect_axis(ds, 'z', where='coords')
    for d in ['time', 'x', 'y', 'z']:
        d = detect_axis(ds, d, where='dims')

    old_depths = srcdomain[depth].values
    new_depths = dstdomain[depth].values

    if depth not in ds.coords:
        ds = ds.assign_coords({depth: ("z", old_depths)})

    output_vars = {}
    for varname in ds.data_vars:
        var = ds[varname]

        if depth not in var.dims:
            continue  # Skip variables without vertical dimension

        if time:
            data = []
            for t in range(var.sizes[time]):
                slice_t = var.isel({time: t}).values
                interp_slice = interpolate(slice_t, old_depths, new_depths, axis=0)
                data.append(interp_slice)
            new_array = np.stack(data, axis=0)
            new_dims = (time, depth) + tuple(d for d in var.dims if d not in [time, depth])
            new_coords = {time: ds[time], "y": ds["y"], "x": ds["x"], depth: new_depths}
            for dim in new_dims:
                if dim not in new_coords:
                    new_coords[dim] = ds[dim]
        else:
            data = var.values
            new_array = interpolate(data, old_depths, new_depths, axis=0)
            new_dims = (depth,) + tuple(d for d in var.dims if d != depth)
            new_coords = {depth: new_depths}
            for dim in new_dims:
                if dim not in new_coords:
                    new_coords[dim] = ds[dim]

        output_vars[varname] = xr.DataArray(new_array, dims=new_dims, coords=new_coords, name=varname, attrs=var.attrs)

    new_ds = xr.merge(output_vars.values())
    new_ds.attrs = ds.attrs
    new_ds.to_netcdf(output_nc)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Interpolate vertical levels using two domain_cfg.nc files.")
    parser.add_argument("--input", "-i", required=True, help="Input NetCDF file to interpolate")
    parser.add_argument("--srcdomain", "-s", required=True, help="Source domain_cfg.nc with original vertical levels")
    parser.add_argument("--dstdomain", "-d", required=True, help="Destination domain_cfg.nc with target vertical levels")
    parser.add_argument("--output", "-o", required=True, help="Output NetCDF file")

    args = parser.parse_args()

    main(args.input, args.srcdomain, args.dstdomain, args.output)

