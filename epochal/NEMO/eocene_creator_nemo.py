
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to create boundary conditions for NEMO EOCENE

"""

import os
import numpy as np
import xarray as xr
from cdo import Cdo

cdo = Cdo()

def make_eocene_tidal_file(herold_folder, output_folder, hcri_value=500.0, cleanup=True):

    m2_varname="eo_tidal_dissipation"

    herold_file_orig = os.path.join(herold_folder, "Green_Huber_eocene_tidal_dissipation_1x1.nc")
    herold_file_remap = os.path.join(herold_folder, "Green_Huber_eocene_tidal_dissipation_r720x360.nc")
    zdfiwm_file = os.path.join(output_folder, "zdfiwm_forcing_r720x360.nc")
    output_file = os.path.join(output_folder, "zdfiwm_forcing_r720x360_eocene.nc")

    if cleanup:
        os.remove(output_file)
        print(f"Cleaning file: {output_file}")

    cdo.remapcon("r720x360", input=herold_file_orig, output=herold_file_remap)

    ds_m2 = xr.open_dataset(herold_file_remap)
    ds_pd = xr.open_dataset(zdfiwm_file)

    m2 = ds_m2[m2_varname]

    # make a copy of present-day file
    ds_out = ds_pd.copy(deep=True)

    # inject M2 into power_cri (1/3 only)
    ds_out['power_cri'][:] = (1.0 / 3.0) * m2.values

    # zero in other energy reservoirs
    ds_out['power_nsq'][:] = 0.0
    ds_out['power_sho'][:] = 0.0
    ds_out['power_bot'][:] = 0.0

    # set vertical length scales
    ds_out['scale_cri'][:] = hcri_value
    ds_out['scale_bot'][:] = 0.0   # unused if power_bot = 0

    # convert NaN into zeroes
    ds_out = ds_out.fillna(0)

    # write output file
    ds_out.to_netcdf(output_file)
    print(f"Written NEMO 4.2 tidal file: {output_file}")

    return None