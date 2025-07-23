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
from cdo import Cdo

# Initialize CDO
cdo = Cdo()


def read_topo(filename):
    """
    Read topography data from a specified path and create land-sea mask, bathymetry, and orography.

    Args:
        path (str): Path to the folder containing the topography files.
    """

    # Check if the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File {filename} not found.")

    # Read the topography data
    data = xr.open_dataset(filename)

    return data


def create_new_topo(path, flag):
    """
    Create a new topographic variable (land-sea mask, opensea mask, bathymetry, orography)
    from the topography data.

    Args:
        data (xarray.Dataset): Topography data.
        path (str): Path to save the output file.
        flag (str): Type of variable to create. One of:
                    "mask_landsea", "mask_opensea", "bathymetry", "orography".
    """

    # Read the topography data
    ds = read_topo(args.path)

    if flag == "mask_landsea":
        ds["landsea_mask"] = (ds["topo"] > 0).astype(int)
        filename = "landsea_mask.nc"

    elif flag == "mask_opensea":
        ds["mask_opensea"] = (ds["topo"] < 0).astype(int)
        filename = "mask_opensea.nc"

    elif flag == "bathymetry":
        ds["bathymetry"] = -ds["topo"].where(ds["topo"] < 0, 0)
        filename = "bathymetry.nc"

    elif flag == "orography":
        ds["orography"] = ds["topo"].where(ds["topo"] > 0, 0)
        filename = "orography.nc"

    else:
        raise ValueError(f"Unknown flag: {flag}")


    # Delete 'topo' from the dataset
    if "topo" in ds:
        del ds["topo"]
    
    # Save to NetCDF
    ds.to_netcdf(os.path.join(path, filename))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Manipulate topography files. ")
    parser.add_argument("-p", "--path", type=str, required=True, help="Path to the topography files.")

    args = parser.parse_args()

    # Extract folder name and file name
    folder_name = os.path.dirname(args.path)
    file_name = os.path.basename(args.path)
    print(f"Folder name: {folder_name}, File name: {file_name}")

    # Create new data variables and save them
    for flag in ["mask_landsea", "mask_opensea", "bathymetry", "orography"]:
        create_new_topo(folder_name, flag)
        print(f"New topography data saved to: {os.path.join(folder_name, f'{flag}.nc')}")
    
    
    cdo.remapcon("N32", input=os.path.join(folder_name, f"landsea_mask.nc"), output=os.path.join(folder_name, f"landsea_mask_remap.nc"))
    cdo.remapcon("N32", input=os.path.join(folder_name, f"orography.nc"), output=os.path.join(folder_name, f"orography_remap.nc"))
    print("All files remapped to N32 grid.")