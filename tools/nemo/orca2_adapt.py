#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tool to modify old ORCA2 file to make them compliant to v4.2.1

Authors
Paolo Davini (CNR-ISAC, Nov 2023)
"""

import xarray as xr

def orca2_fixer(data):

    """Inner fixer function to be applied to each dataarray"""

# Check if the x and y dimensions exist
    if 'x' in data.dims and 'y' in data.dims:

        # Clean borders
        data = data.isel(x=slice(1, -1), y=slice(None, -1))

        print(f"Processing completed for {data.name}")
    else:
        print(f"{data.name} does not have the required x and y dimensions.")

    return data

def orca2_main(input, output):

    """Main fixer function to be applied to each dataset"""

    ds = xr.open_dataset(input, decode_times=False)
    newds = ds.map(orca2_fixer)
    encoding = {var: {'_FillValue': None} for var in newds.data_vars}
    newds.to_netcdf(output, encoding=encoding)

    return newds


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.nc output_file.nc")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    orca2_main(input_file, output_file)

