#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A tool to produce the mask_util from mesh_mask
Also modifies the domain_cfg file by renaming the vertical dimension

Authors
Paolo Davini (CNR-ISAC, Jul 2025)
"""

import argparse
import os
import xarray as xr


mesh_file = 'mesh_mask.nc'
domain_file = 'domain_cfg.nc'

default_src_dir = '/lus/h2resw01/hpcperm/ccpd/ecearth4/revisions/main/sources/nemo-4.2/tools/DOMAINcfg'
default_tgt_dir = '/lus/h2resw01/hpcperm/ccpd/PALEORCA2'

def domain_cfg(src_dir, tgt_dir):
    """
    Modify the domain_cfg file by renaming the vertical dimension and resetting it. 

    Args:
    src_dir (str): Path to the source directory containing the domain_cfg file. 
    tgt_dir (str): Path to the target directory where the modified domain_cfg file will be saved.
    
    """

    # load the xarray files
    domain = xr.open_dataset(f'{src_dir}/{domain_file}')

    # rename and reset the vertical dimension
    domain = domain.rename_dims({'nav_lev': 'z'})
    domain = domain.reset_index('nav_lev').reset_coords('nav_lev')


    # set the fill values
    encoding_var = {var: {'_FillValue': None} for var in domain.data_vars}
    encoding_coord = {var: {'_FillValue': None} for var in domain.coords}
    encoding = {**encoding_var, **encoding_coord}

    # write the file
    os.makedirs(tgt_dir, exist_ok=True)
    domain.to_netcdf(f'{tgt_dir}/domain_cfg.nc', encoding=encoding, unlimited_dims=['time_counter'])

def maskutil(src_dir, tgt_dir):
    """
    Extracts mask variables from a mesh dataset and saves them to a new maskutil.nc file.

    Args:
    src_dir (str): Path to the source directory containing the mesh_mask.nc file.
    tgt_dir (str): Path to the target directory where the maskutil.nc file will be saved.

    Returns:
    None
    """
    mesh = xr.open_dataset(f'{src_dir}/mesh_mask.nc')
    masks = mesh[['tmaskutil','umaskutil','vmaskutil']]
    masks = masks.rename_dims({'time_counter': 't'}).drop_vars('time_counter')
    masks.attrs = {'Conventions': "CF-1.1"}
    os.makedirs(tgt_dir, exist_ok=True)
    masks.to_netcdf(f'{tgt_dir}/maskutil.nc', unlimited_dims=['t'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A tool to modify domain_cfg ORCA2 file')



    parser.add_argument('--src_dir', type=str, default=default_src_dir, 
                        help='Path to src domain directory (domainCFG path)')
    parser.add_argument('--tgt_dir', type=str, default=default_tgt_dir, 
                        help='Path to target directory')

    args = parser.parse_args()

    domain_cfg(args.src_dir, args.tgt_dir)
    maskutil(args.src_dir, args.tgt_dir)
