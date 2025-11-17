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

def vegetation_zhang(self, field):
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

        tvh = xr.full_like(field['tvh'], 0)
        tvl = xr.full_like(field['tvl'], 0)
        cvh = xr.zeros_like(tvh)
        cvl = xr.zeros_like(tvl)

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

            if biome_id in biome_to_tvh:
                tvh = xr.where(mask, biome_to_tvh[biome_id], tvh)
                cvh = xr.where(mask, 1.0, cvh)
            elif biome_id in biome_to_tvl:
                tvl = xr.where(mask, biome_to_tvl[biome_id], tvl)
                cvl = xr.where(mask, 1.0, cvl)
            else:
                print(f"Warning: biome {biome_id} not in mapping.")

        # Replace field variables
        for vname, newval in zip(["tvh", "tvl", "cvh", "cvl"], [tvh, tvl, cvh, cvl]):
            field[vname].data = newval.data

        return field

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
    

def aerosols(self):
    """
    Convert Herold Eocene aerosol data (kg/kg) to column-integrated mass (kg/m²)
    on the IFS grid, using the US Standard Atmosphere for density interpolation.
    """

    print("→ Loading Herold aerosol data")
    paleoaerfile = os.path.join(self.herold, "herold_etal_eocene_CAM4_BAM_aerosols.nc")
    aer_paleo = xr.open_dataset(paleoaerfile)

    # Load IFS reference aerosol climatology
    aerfile_ifs_pd = os.path.join(self.idir, 'oifs', 'ifsdata/aerosol_cams_climatology_43R3a.nc')
    aer_ifs_pd = xr.open_dataset(aerfile_ifs_pd)
    aer_ifs_paleo = aer_ifs_pd.copy()

    # Re-bin the Herold aerosol bins
    aer_paleo_newbin = aer_paleo.copy()
    aer_paleo_newbin['DST01'] = aer_paleo['DST01']*0.57
    aer_paleo_newbin['DST02'] = aer_paleo['DST01']*0.39
    aer_paleo_newbin['DST03'] = aer_paleo['DST02'] + aer_paleo['DST03'] + aer_paleo['DST04']
    aer_paleo_newbin['SSLT01'] = aer_paleo['SSLT01']*0.57
    aer_paleo_newbin['SSLT02'] = aer_paleo['SSLT01']*0.39 + aer_paleo['SSLT02'] + aer_paleo['SSLT03']
    aer_paleo_newbin['SSLT03'] = aer_paleo['SSLT04']

    # Check or generate US Standard Atmosphere
    ua_file = os.path.join(self.herold, "us_standard_atmosphere_newlevs.nc")
    if not os.path.exists(ua_file):
        print("→ Generating US Standard Atmosphere data")
        # US Standard Atmosphere data (alt, temp, g, pressure, density, mu)
        data = """
        -1000   21.50   9.810   11.39   1.347   1.821
        0       15.00   9.807   10.13   1.225   1.789
        1000    8.50    9.804   8.988   1.112   1.758
        2000    2.00    9.801   7.950   1.007   1.726
        3000    -4.49   9.797   7.012   0.9093  1.694
        4000    -10.98  9.794   6.166   0.8194  1.661
        5000    -17.47  9.791   5.405   0.7364  1.628
        6000    -23.96  9.788   4.722   0.6601  1.595
        7000    -30.45  9.785   4.111   0.5900  1.561
        8000    -36.94  9.782   3.565   0.5258  1.527
        9000    -43.42  9.779   3.080   0.4671  1.493
        10000   -49.90  9.776   2.650   0.4135  1.458
        15000   -56.50  9.761   1.211   0.1948  1.422
        20000   -56.50  9.745   0.5529  0.08891 1.422
        25000   -51.60  9.730   0.2549  0.04008 1.448
        30000   -46.64  9.715   0.1197  0.01841 1.475
        40000   -22.80  9.684   0.0287  0.003996 1.601
        50000   -2.5    9.654   0.007978 0.001027 1.704
        60000   -26.13  9.624   0.002196 0.0003097 1.584
        70000   -53.57  9.594   0.00052  0.00008283 1.438
        80000   -74.51  9.564   0.00011  0.00001846 1.321
        """
        lines = data.strip().split('\n')
        parsed_data = [list(map(float, line.split())) for line in lines]
        arr = np.array(parsed_data)
        altitude = arr[:,0]
        temperature = 273.15 + arr[:,1]
        g = arr[:,2]
        pressure = 1e2*arr[:,3]
        density = arr[:,4]
        mu = arr[:,5]

        ds = xr.Dataset(
            {
                'temperature': (['altitude'], temperature),
                'gravity': (['altitude'], g),
                'pressure': (['altitude'], pressure),
                'density': (['altitude'], density),
                'viscosity': (['altitude'], mu)
            },
            coords={'altitude': (['altitude'], altitude)}
        )

        # log-pressure coordinate & interpolate to Herold levels
        log_pressure = np.log(ds.pressure)
        ds = ds.assign_coords(log_pressure=log_pressure)
        ds_new = ds.swap_dims({'altitude':'log_pressure'})
        new_log_pressure = np.log(aer_paleo.lev)
        ds_interp = ds_new.interp(log_pressure=new_log_pressure, method='linear')
        ds_interp.to_netcdf(ua_file)
        print(f"→ US Standard Atmosphere saved at {ua_file}")
    else:
        ds_interp = xr.open_dataset(ua_file)

    # Map Herold -> IFS variable names
    var_dict = {
        'Sulfates': "SO4",
        'Black_Carbon_hydrophilic': "CB1",
        'Black_Carbon_hydrophobic': "CB2",
        'Mineral_Dust_bin1': "DST01",
        'Mineral_Dust_bin2': "DST02",
        'Mineral_Dust_bin3': "DST03",
        'Organic_Matter_hydrophilic': "OC1",
        'Organic_Matter_hydrophobic': "OC2",
        'Sea_Salt_bin1': "SSLT01",
        'Sea_Salt_bin2': "SSLT02",
        'Sea_Salt_bin3': "SSLT03"
    }

    # Convert to column mass and regrid
    for varname2 in var_dict:
        varname1 = var_dict[varname2]
        print(f"→ Processing {varname1} -> {varname2}")
        new_var = aer_paleo_newbin[varname1] * ds_interp.density
        new_var_rg = regrid_dataset(new_var, regrid_to_reference=aer_ifs_paleo[varname2])
        new_var_rg_int = -new_var_rg.integrate(coord='altitude')
        aer_ifs_paleo[varname2].data = new_var_rg_int.data.astype('float32')

    # Save Eocene aerosol climatology
    output_path = os.path.join(self.odir_init, "aerosol_cams_climatology_43R3a_EOCENO.nc")
    if os.path.exists(output_path):
        os.remove(output_path)
    aer_ifs_paleo.to_netcdf(output_path)
    print(f"→ Eocene aerosol data saved at {output_path}")
    return output_path

def compute_slope (grib_field, sd_eoc, a=4.376786e-05, b=2.476405e-04):
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
    # Interpolate sd_eoc to the GRIB grid if needed
    sd = sd_eoc

    # Apply the linear function
    new_slope = a * sd + b

    return new_slope

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