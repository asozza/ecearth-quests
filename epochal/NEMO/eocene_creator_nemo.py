
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to create boundary conditions for NEMO EOCENE

"""

import os
import numpy as np
import xarray as xr
import shutil
import tempfile
import subprocess
import pandas as pd
from cdo import Cdo

cdo = Cdo()

class EoceneNEMO():

    def __init__(self, herold_folder, input_folder, output_folder, grid= "r720x360"):
        self.herold_folder = herold_folder
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.cdo = Cdo()

        self.grid=grid

    def create_tidal_mixing(self, hcri_value=500.0, cleanup=True):

        m2_varname="eo_tidal_dissipation"

        herold_file_orig = os.path.join(self.herold_folder, "Green_Huber_eocene_tidal_dissipation_1x1.nc")
        herold_file_remap = os.path.join(self.herold_folder, "Green_Huber_eocene_tidal_dissipation_r720x360.nc")
        zdfiwm_file = os.path.join(self.input_folder, "nemo/initial/zdfiwm_forcing_r720x360.nc")
        output_file = os.path.join(self.output_folder, "nemo/initial/zdfiwm_forcing_r720x360.nc")

        if cleanup and os.path.exists(output_file):
            os.remove(output_file)
            print(f"Cleaning file: {output_file}")

        self.cdo.remapcon(self.grid, input=herold_file_orig, output=herold_file_remap)

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

        ds_m2.close()
        ds_pd.close()

        return None
    
    def create_geothermal_flux(self, cap_value=400, cleanup=True):

        gh_file = os.path.join(self.input_folder,"nemo/initial/Goutorbe_ghflux.nc")
        age_file = os.path.join(self.herold_folder,"age_depth_bath_55.xyadb")
        gh_file_out = os.path.join(self.output_folder,"nemo/initial/Goutorbe_ghflux.nc")

        if cleanup and os.path.exists(gh_file_out):
            os.remove(gh_file_out)
            print(f"Cleaning file: {gh_file_out}")
        
        # Loading Geothermal flux NEMO file
        ds_gh = xr.open_dataset(gh_file)
        lat_tgt = ds_gh.lat.values
        lon_tgt = ds_gh.lon.values

        age = pd.read_csv(age_file, sep=r"\s+", header=None, usecols=[0,1,2], names=["lon","lat","age"], engine="c")
        t = age.set_index(["lat", "lon"])["age"].to_xarray()

        q = xr.where(
            t <= 55,
            510 * t**(-0.5),
            48 + 96 * np.exp(-0.0278 * t)
        )
        
        q_min = float(q.min())
        q = q.fillna(q_min)
    
        # Flip lat if needed (so it’s increasing)
        q = q.sortby("lat")

        # 1D interpolation along lon
        tmp = np.empty((q.lat.size, lon_tgt.size))
        for i in range(q.lat.size):
            tmp[i, :] = np.interp(lon_tgt, q.lon.values, q.values[i, :])

        # 1D interpolation along lat
        q_interp_values = np.empty((lat_tgt.size, lon_tgt.size))
        for j in range(len(lon_tgt)):
            q_interp_values[:, j] = np.interp(lat_tgt, q.lat.values, tmp[:, j])

        # Wrap as DataArray on NEMO grid
        q_interp = xr.DataArray(
            q_interp_values,
            coords={"lat": lat_tgt, "lon": lon_tgt},
            dims=("lat", "lon"),
            name="gh_flux"
        )

        # Cap values
        q_interp = xr.where(q_interp > cap_value, cap_value, q_interp)

        if lat_tgt[0] > lat_tgt[-1]:
            q_interp = q_interp.sortby("lat", ascending=False)
        # Assign into Goutorbe_ghflux.nc
        ds_gh["gh_flux"][0, :, :] = q_interp

        
        ds_gh.to_netcdf(gh_file_out)
        print(f"Written NEMO 4.2 tidal file: {gh_file_out}")
        ds_gh.close()

        return None