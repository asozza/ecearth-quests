
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to create boundary conditions for NEMO EOCENE

"""

import os
import xarray as xr
import numpy as np
import shutil
from cdo import Cdo

cdo = Cdo()

import os
import xarray as xr
from cdo import Cdo

cdo = Cdo()


def make_eocene_tidal_file(herold_folder, output_folder, hcri_value=500.0):
    """ 
    Create IWM forcing file compliant with NEMO 4.2, according to de Lavergne setup:
        - 4 fields: power_cri, power_nsq, power_sho, power_bot
        - 2 length scales: scale_cri, scale_bot
    
    Usage:

    make_eocene_tidal_file(
        herold_folder="/perm/itas/data/deepMIP/Herold_etal_2014", 
        output_folder="/ec/res4/hpcperm/itas/data/ece-4-database/nemo/initial", 
        hcri_value=500.0
    )
    """

    m2_varname="eo_tidal_dissipation"

    herold_file_orig = os.path.join(herold_folder, "Green_Huber_eocene_tidal_dissipation_1x1.nc")
    herold_file_remap = os.path.join(herold_folder, "Green_Huber_eocene_tidal_dissipation_r720x360.nc")
    zdfiwm_file = os.path.join(output_folder, "zdfiwm_forcing_r720x360.nc")
    output_file = os.path.join(output_folder, "zdfiwm_forcing_r720x360_eocene.nc")

    cdo.remapcon("r720x360", input=herold_file_orig, output=herold_file_remap)

    ds_m2 = xr.open_dataset(herold_file_remap)
    ds_pd = xr.open_dataset(zdfiwm_file)

    m2 = ds_m2[m2_varname]

    # --- deep copy of PD file
    ds_out = ds_pd.copy(deep=True)

    # --- inject M2 into power_cri (1/3 only)
    ds_out['power_cri'][:] = (1.0 / 3.0) * m2.values

    # --- zero other energy reservoirs
    ds_out['power_nsq'][:] = 0.0
    ds_out['power_sho'][:] = 0.0
    ds_out['power_bot'][:] = 0.0

    # --- vertical length scales
    ds_out['scale_cri'][:] = hcri_value
    ds_out['scale_bot'][:] = 0.0   # unused if power_bot = 0

    # --- write file
    ds_out.to_netcdf(output_file)

    print(f"Written NEMO 4.2 tidal file: {output_file}")
