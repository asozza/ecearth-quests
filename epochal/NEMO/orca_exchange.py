#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to exchange fields between NEMO datasets, as along as they share the same dimensions and structure.

Author: Alessandro Sozza (CNR-ISAC)
Date: July 2025
"""

# example:
# ./orca_levels.py -i /ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial/woa13-levitus.nc -s /ec/res4/hpcperm/itas/data/ece-4-database/nemo/domain/eORCA1/domain_cfg.nc -d /ec/res4/hpcperm/itas/data/ece-4-database/nemo/domain/ORCA2/domain_cfg.nc -o /ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial/woa13-levitus_ORCA2.nc


import sys
import os
import argparse
import xarray as xr
import numpy as np
from scipy.interpolate import interp1d

def main(input_nc, output_nc):
    
    ds_in= xr.open_dataset(input_nc)
    ds_out = xr.open_dataset(output_nc)

    for var in ['lat', 'lon']
    
    
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

