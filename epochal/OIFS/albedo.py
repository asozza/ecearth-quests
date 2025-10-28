"""Some utilities for OIFS grid definition"""
import re
import os
import tempfile
import shutil
import numpy as np
import xarray as xr
import subprocess
from cdo import Cdo
cdo = Cdo()

def albedo(field: xr.Dataset, var=None, lsm_present=None, landsea=None, **kwargs):
    """
    Apply both:
    Land-sea mask-based albedo reconstruction using `lsm_present`
    Eocene land-sea mask adjustment (albedo=0.05, LAI=NaN over ocean)
    """

    # VALIDATION
    if lsm_present is None:
        raise ValueError("You must provide `lsm_present` (present-day land-sea mask).")
    if landsea is None:
        raise ValueError("You must provide `landsea` (Eocene land-sea mask).")

    print("Applying combined Albedo reconstruction + Eocene mask...")

    # PREPROCESS MASKS
    # Ensure ascending latitude
    if not np.all(np.diff(lsm_present['lat']) > 0):
        lsm_present = lsm_present.sortby('lat')
    if not np.all(np.diff(landsea['lat']) > 0):
        landsea = landsea.sortby('lat')

    # Extract variable if landsea is a dataset
    if isinstance(landsea, xr.Dataset):
        mask_var = list(landsea.data_vars)[0]
        landsea = landsea[mask_var]

    # Interpolate Eocene mask to field grid
    landsea_interp = landsea.interp(lat=field["lat"], lon=field["lon"], method="nearest")

    # Match time dimension if needed
    if "time" in field.dims:
        ntime = field.dims["time"]
        if "time" not in landsea_interp.dims or landsea_interp.sizes.get("time", 1) != ntime:
            landsea_interp = landsea_interp.expand_dims("time").broadcast_like(field.isel(time=slice(0, 1)))
            landsea_interp = xr.concat([landsea_interp] * ntime, dim="time")

    # Boolean masks
    present_mask = (lsm_present.broadcast_like(field)) > 0.5
    eocene_mask = (landsea_interp > 0.5).transpose(*field[list(field.data_vars)[0]].dims)

    # APPLY ALBEDO RECONSTRUCTION
    for v in field.data_vars:
        da = field[v]

        # Apply the present-day mask (land only)
        masked = da.where(present_mask)

        # Sort latitude ascending
        flip = False
        if da["lat"].values[0] > da["lat"].values[-1]:
            da = da.sortby("lat")
            masked = masked.sortby("lat")
            flip = True

        # Zonal mean over land
        zonal_mean = masked.mean(dim="lon", skipna=True)
        zonal_mean_filled = zonal_mean.interpolate_na(dim="lat", method="nearest")
        da_recon = zonal_mean_filled.broadcast_like(da)

        # Polar band filling
        lat = da["lat"]
        if (lat < -52).any():
            mean_s = da_recon.sel(lat=lat.where((lat >= -52) & (lat <= -46), drop=True)).mean(dim=("lat", "lon"), skipna=True)
            da_recon = da_recon.where(~(lat < -52), other=mean_s)
        if (lat > 75).any():
            mean_n = da_recon.sel(lat=lat.where((lat >= 70) & (lat <= 75), drop=True)).mean(dim=("lat", "lon"), skipna=True)
            da_recon = da_recon.where(~(lat > 75), other=mean_n)

        if flip:
            da_recon = da_recon.sortby("lat", ascending=False)

        # Replace variable data
        field[v].data = da_recon.data

    print("Albedo reconstruction complete.")

    # APPLY EOCENE MASK RULES
    albedo_vars = ["al", "aluvp", "aluvd", "alnip", "alnid"]
    lai_vars = ["lai_lv", "lai_hv"]

    for v in albedo_vars:
        if v in field:
            field[v].data = np.where(eocene_mask, field[v].data, 0.05)

    for v in lai_vars:
        if v in field:
            field[v].data = np.where(eocene_mask, field[v].data, np.nan)

    print("Eocene land-sea mask applied successfully.")
    print("Combined modification complete, GRIB structure preserved.")
    return field