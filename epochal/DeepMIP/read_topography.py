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


def create_new_topo(data, path, flag):
    """
    Create a new topographic variable (land-sea mask, opensea mask, bathymetry, orography)
    from the topography data.

    Args:
        data (xarray.Dataset): Topography data.
        path (str): Path to save the output file.
        flag (str): Type of variable to create. One of:
                    "mask_landsea", "mask_opensea", "bathymetry", "orography".
    """

    if flag == "mask_landsea":
        data["landsea_mask"] = (data["topo"] > 0).astype(int)
        filename = "landsea_mask.nc"

    elif flag == "mask_opensea":
        data["mask_opensea"] = (data["topo"] < 0).astype(int)
        filename = "mask_opensea.nc"

    elif flag == "bathymetry":
        data["bathymetry"] = -data["topo"].where(data["topo"] < 0)
        filename = "bathymetry.nc"

    elif flag == "orography":
        data["orography"] = data["topo"].where(data["topo"] > 0)
        filename = "orography.nc"

    else:
        raise ValueError(f"Unknown flag: {flag}")

    # Save to NetCDF
    data.to_netcdf(os.path.join(path, filename))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Manipulate topography files. ")
    parser.add_argument("-p", "--path", type=str, required=True, help="Path to the folder containing the topography files.")

    args = parser.parse_args()

    # Read the topography data
    data = read_topo(args.path)

    # Create new data variables and save them
    for flag in ["mask_landsea", "mask_opensea", "bathymetry", "orography"]:
        create_new_topo(data, args.path, flag)
        print(f"New topography data saved to: {os.path.join(args.path, f'{flag}.nc')}")
    