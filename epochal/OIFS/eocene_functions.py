"""
Functions to compute Eocene-specific fields:
- slope from SD
- remapped Herold fields
- vegetation mapping
- albedo transformations
- aerosol reconstructions
"""

import re
import os
import tempfile
import shutil
import numpy as np
import xarray as xr
import subprocess
import xesmf as xe
import shutil
import tempfile
from cdo import Cdo
cdo = Cdo()

def _to_tll(field):
    """Force (time, lat, lon) ordering"""
    if set(field.dims) == {"time", "lat", "lon"}:
        return field.transpose("time", "lat", "lon")
    return field

def vegetation_zhang(field, var=None, herold_path=None, gaussian=None, **kwargs):
        """"
        Alternative method to create the ICMGG vegetation data for the Eocene OIFS.
        Replace the vegetation data with the one from the Herold data.
        Set the vegetation content to 1 for the dominant vegetation type and 0 to the others.
        Perform a mapping using the Zhang et al., 2021 criteria. 
        Always returns a Dataset with tvh, tvl, cvh, cvl.
        """
        # --- Load Herold biome data and remap ---
        herold_file = os.path.join(herold_path, "herold_etal_eocene_biome_1x1.nc")
        herold_remap = cdo.remapnn(
            f"N{gaussian}", input=herold_file, 
            output=os.path.join(herold_path, "herold_etal_eocene_biome_1x1_N32.nc")
            )

        herold = xr.open_dataset(herold_remap)

        tvh = _to_tll(xr.full_like(field["tvh"], 0))
        tvl = _to_tll(xr.full_like(field["tvl"], 0))

        cvh = xr.zeros_like(tvh)
        cvl = xr.zeros_like(tvl)

        print("Initial shapes:")
        print(f"tvh: {tvh.shape}, tvl: {tvl.shape}, cvh: {cvh.shape}, cvl: {cvl.shape}")

        # === Biome to vegetation ID mappings ===
        biome_to_tvh = {
            1: 1, # Tropical forest → Evergreen broadleaf trees
            2: 2, # Warm-temperate forest → Evergreen needleleaf trees
            6: 3, # Temperate forest → Deciduous broadleaf
            7: 4, # Boreal forest → Deciduous needleleaf
        }

        biome_to_tvl = {
            3: 5, # Savanna → Tall grass
            4: 6, # Grassland → Short grass
            5: 7, # Desert → Semidesert
            8: 8, # Tundra → Tundra
            9: 8, # Dry Tundra → Tundra
        }

        for biome_id in range(1, 10): # assuming biome IDs go from 1 to 9
            mask = herold['eocene_biome_hp'] == biome_id
            # Expand mask to time dimension if present
            if 'time' in tvh.dims:
                mask = mask.expand_dims(time=tvh['time'])
            # Sanity check
            print(f"Biome {biome_id}: mask shape {mask.shape}")
            for arr_name, arr in zip(['tvh','tvl'], [tvh, tvl]):
                print(f"Array {arr_name} shape: {arr.shape}")
                if mask.shape != arr.shape:
                    print(f"WARNING: mask shape {mask.shape} != {arr_name} shape {arr.shape}")
    

            if biome_id in biome_to_tvh:
                tvh = xr.where(mask, biome_to_tvh[biome_id], tvh)
                cvh = xr.where(mask, 1.0, cvh)
            elif biome_id in biome_to_tvl:
                tvl = xr.where(mask, biome_to_tvl[biome_id], tvl)
                cvl= xr.where(mask, 1.0, cvl)
            else:
                print(f"Warning: biome {biome_id} not in mapping.")

        field = xr.Dataset({
            "tvh": tvh,
            "tvl": tvl,
            "cvh": cvh,
            "cvl": cvl,
        }, coords=field.coords, attrs=field.attrs)
        
        print(field)
        return field

        # Sanity check
        #print(f"Biome {biome_id}: mask shape {mask.shape}")
        #for arr_name, arr in zip(['tvh','tvl'], [tvh, tvl]):
        #    print(f"Array {arr_name} shape: {arr.shape}")
         #   if mask.shape != arr.shape:
         #       print(f"WARNING: mask shape {mask.shape} != {arr_name} shape {arr.shape}")
    

        #    if biome_id in biome_to_tvh:
        #        tvh = xr.where(mask, biome_to_tvh[biome_id], tvh)
        #        cvh = xr.where(mask, 1.0, cvh)
        #    elif biome_id in biome_to_tvl:
        #        tvl = xr.where(mask, biome_to_tvl[biome_id], tvl)
        #        cvl = xr.where(mask, 1.0, cvl)
        #    else:
        #        print(f"Warning: biome {biome_id} not in mapping.")

        # Final shape check before assignment
        #for var_name, arr in zip(["tvh","tvl","cvh","cvl"], [tvh, tvl, cvh, cvl]):
        #    print(f"Final {var_name} shape: {arr.shape}, field shape: {field[var_name].shape}")
        #    if arr.shape != field[var_name].shape:
        #        print(f"WARNING: {var_name} shape mismatch! {arr.shape} vs {field[var_name].shape}")

        # Assign back to field
        for name, arr in data_vars.items():
            field[name].data = arr.data

        return field
        # Replace field variables
        #for var, newval in zip(["tvh", "tvl", "cvh", "cvl"], [tvh, tvl, cvh, cvl]):
        #    field[var].data = newval.data

        #return field

def prepare_vegetation_zhang(self):
        """"
        Alternative method to create the ICMGG vegetation data for the Eocene OIFS.
        Replace the vegetation data with the one from the Herold data.
        Set the vegetation content to 1 for the dominant vegetation type and 0 to the others.
        Perform a mapping using the Zhang et al., 2021 criteria. 
        """


        herold_file = os.path.join(self.herold, "herold_etal_eocene_biome_1x1.nc")
        herold_remap = cdo.remapnn(
            f"N{self.gaussian}", 
            input=herold_file, 
            output=os.path.join(self.herold, "herold_etal_eocene_biome_1x1_N32.nc")
        )

        herold = xr.open_dataset(herold_remap)

        # Initialize arrays with the correct shape and coordinates
        tvh = _to_tll(xr.DataArray(np.zeros_like(field['tvh'].values), dims=field['tvh'].dims, coords=coords))
        tvl = _to_tll(xr.DataArray(np.zeros_like(field['tvl'].values), dims=field['tvl'].dims, coords=coords))
        cvh = xr.DataArray(np.zeros_like(tvh.values), dims=tvh.dims, coords=coords)
        cvl = xr.DataArray(np.zeros_like(tvl.values), dims=tvl.dims, coords=coords)

        # === Biome to vegetation ID mappings ===
        biome_to_tvh = {
            1: 1, # Tropical forest → Evergreen broadleaf trees
            2: 2, # Warm-temperate forest → Evergreen needleleaf trees
            6: 3, # Temperate forest → Deciduous broadleaf
            7: 4, # Boreal forest → Deciduous needleleaf
        }

        biome_to_tvl = {
            3: 5, # Savanna → Tall grass
            4: 6, # Grassland → Short grass
            5: 7, # Desert → Semidesert
            8: 8, # Tundra → Tundra
            9: 8, # Dry Tundra → Tundra
        }

        # === Create blank data arrays ===
        shape = herold['eocene_biome_hp'].shape
        coords = herold.coords

        tvh = xr.full_like(herold['eocene_biome_hp'], fill_value=0)
        tvl = xr.full_like(herold['eocene_biome_hp'], fill_value=0)
        cvh = xr.zeros_like(tvh)
        cvl = xr.zeros_like(tvl)

        for biome_id in range(1, 10): # assuming biome IDs go from 1 to 9
            mask = herold['eocene_biome_hp'] == biome_id

            if biome_id in biome_to_tvh:
                tvh = xr.where(mask, biome_to_tvh[biome_id], tvh)
                cvh = xr.where(mask, 1.0, cvh)
            elif biome_id in biome_to_tvl:
                tvl = xr.where(mask, biome_to_tvl[biome_id], tvl)
                cvl = xr.where(mask, 1.0, cvl)
            else:
                print(f"Warning: biome {biome_id} not in mapping.")

        # === Assemble final dataset ===
        vegetation_ds = xr.Dataset(
        {
            "tvh": tvh,
            "tvl": tvl,
            "cvh": cvh,
            "cvl": cvl
        },

        )

        # === Save ===
        output_path = os.path.join(self.odir_init, "ICMGG_vegetation.nc")
        if os.path.exists(output_path):
            os.remove(output_path)
        vegetation_ds.to_netcdf(output_path)

        return output_path


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
            field[v].data = np.where(eocene_mask, field[v].data, 0)

    print("Eocene land-sea mask applied successfully.")
    print("Combined modification complete, GRIB structure preserved.")
    return field

def compute_slope (field, var=None, sd_eoc=None, a=4.376786e-05, b=2.476405e-04):
    """
    Replace a GRIB field (slor) using the linear transfer function:
        slope = a * sd + b
    Parameters
    ----------
    grib_field : xarray.DataArray
        The existing slope field from GRIB (ignored except for shape alignment).
    sd_eoc : xarray.DataArray
        Eocene sd_orography on the model grid.
    a, b : float
        Linear transfer coefficients.
    """
    
    # Compute slope
    new_slope = a * sd_eoc + b
    print("Computed new slope (before expanding dims):", new_slope.shape)

    for v in var:
        target = field[v]

        # Ensure time dimension exists and shape matches
        if 'time' in target.dims and len(target.dims) == 3:
            # GRIB has time dimension (1, lat, lon)
            if new_slope.shape != target.shape:
                # expand time dim if missing
                new_slope_exp = new_slope
                if 'time' not in new_slope_exp.dims:
                    new_slope_exp = new_slope_exp.expand_dims('time', axis=0)
                # broadcast to match exactly
                if new_slope_exp.shape != target.shape:
                    new_slope_exp = np.broadcast_to(new_slope_exp, target.shape)
                new_slope_to_assign = new_slope_exp
            else:
                new_slope_to_assign = new_slope
        else:
            new_slope_to_assign = new_slope

        print(f"Assigning slope to variable '{v}' with shape {target.shape}")
        field[v].data = new_slope_to_assign.data

    return field

def prepare_vegetation(self):
    """"
    Create the ICMGG vegetation data for the Eocene OIFS.
    Replace the vegetation data with the one from the Herold data.
    Set the vegetation type to 0 for all types.
    Perform a mapping from present-day initial conditions
    """


    herold_file = os.path.join(self.herold, "herold_etal_eocene_biome_1x1.nc")
    herold_remap = cdo.remapnn(
        f"N{self.gaussian}", 
        input=herold_file, 
        output=os.path.join(self.herold, "herold_etal_eocene_biome_1x1_N32.nc")
    )
            
    icmgg_file = os.path.join(self.idir_init, "ICMGGECE4INIT")
    if os.path.exists(os.path.join(self.herold, "ICMGG.nc")):
        os.remove(os.path.join(self.herold, "ICMGG.nc"))
    icmgg_remap = cdo.setgridtype(
        "regularnn", 
        input=icmgg_file, 
        output=os.path.join(self.herold ,"ICMGG.nc"),
        options=NC4
    )

    herold = xr.open_dataset(herold_remap)
    icmgg = xr.open_dataset(icmgg_remap)

    biome_dict = {'tvh': {}, 'tvl': {}}
    for vegtype in ["tvh", "tvl"]:
        for i in range(1, 11):
            vegid = icmgg[vegtype].where(herold["prei_biome_hp"] == i).values
            vegid = vegid[~np.isnan(vegid)]
            unique, counts = np.unique(vegid, return_counts=True)
            if unique.size>0:
                biome_dict[vegtype][i] = int(unique[np.argmax(counts)])
            else:
                biome_dict[vegtype][i] = None

    eocene_icmgg = icmgg[['tvh', 'tvl', 'cvh', 'cvl']]
    for vegtype in ["tvh", "tvl"]:
        eocene_icmgg[vegtype] = eocene_icmgg[vegtype]*0
        for i in range(1, 11):
            eocene_icmgg[vegtype] = xr.where(
                herold['eocene_biome_hp'] == i,
                biome_dict[vegtype][i],
                eocene_icmgg[vegtype])
        
    for vegtype in ["cvh", "cvl"]:
        eocene_icmgg[vegtype] = eocene_icmgg[vegtype]*0

    if os.path.exists(os.path.join(self.odir_init, "ICMGG_vegetation.nc")):
        os.remove(os.path.join(self.odir_init, "ICMGG_vegetation.nc"))
    eocene_icmgg.to_netcdf(
        os.path.join(self.odir_init, "ICMGG_vegetation.nc")
    )

    return os.path.join(self.odir_init, "ICMGG_vegetation.nc")