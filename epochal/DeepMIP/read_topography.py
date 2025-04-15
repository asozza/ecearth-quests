#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to read topography from Herold et al. 2014 and create:
- land-sea mask
- mask_opensea
- bathymetry
- orography

Usage: ./read_topography.py --path <path>

    -h, --help                Show this help message and exit.
    -p, --path <folder-path>  Path to the folder containing the topography files.

Authors
Alessandro Sozza (CNR-ISAC, April 2025)
"""

import os
import argparse
import xarray as xr


def read_topo(path):
    """
    Read topography data from a specified path and create land-sea mask, bathymetry, and orography.

    Args:
        path (str): Path to the folder containing the topography files.
    """

    # Check if the folder exists
    if not os.path.isdir(args.path):
        raise FileNotFoundError(f"The specified folder does not exist: {args.path}")

    # Read the topography data
    data = xr.open_dataset(os.path.join(path, "topo.nc"))

    return data

def create_landsea_mask(data, path):
    """
    Create land-sea mask, bathymetry, and orography from the topography data.
    Args:
        data (xarray.Dataset): Topography data.       
    """

    # Create land-sea mask
    data['landsea_mask'] = (data['topo'] > 0).astype(int)

    # Save the results
    data.to_netcdf(os.path.join(path, "landsea_mask.nc"))

def create_opensea_mask(data, path):
    """
    Create open sea mask from the topography data.
    Args:
        data (xarray.Dataset): Topography data.
    """
    # Create open sea mask
    data['mask_opensea'] = (data['topo'] < 0).astype(int)

    # Save the results
    data.to_netcdf(os.path.join(path, "mask_opensea.nc"))

def create_bathymetry(data, path):
    """
    Create bathymetry from the topography data.
    Args:
        data (xarray.Dataset): Topography data.
    """

    # Create bathymetry
    data['bathymetry'] = -data['topo'].where(data['topo'] < 0)

    # save the results
    data.to_netcdf(os.path.join(path, "bathymetry.nc"))

def create_orography(data, path):
    """
    Create orography from the topography data.
    Args:
        data (xarray.Dataset): Topography data.
    """

    # Create orography
    data['orography'] = data['topo'].where(data['topo'] > 0)

    # save the results
    data.to_netcdf(os.path.join(path, "orography.nc"))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Manipulate topography files. ")
    parser.add_argument("-p", "--path", type=str, required=True, help="Path to the folder containing the topography files.")

    args = parser.parse_args()

    # Read the topography data
    data = read_topo(args.path)

    # Create new data with land-sea mask, bathymetry, and orography
    create_landsea_mask(data, args.path)
    create_opensea_mask(data, args.path)
    create_bathymetry(data, args.path)
    create_orography(data, args.path)
    print(f"New topography data saved to: f{args.path}")
    