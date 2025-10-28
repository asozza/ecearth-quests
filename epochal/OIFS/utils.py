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

# CDO grib2 options
GRIB2="-f grb2 --eccodes"
GRIB1="-f grb1 --eccodes"
NC4='-f nc4 --eccodes'


def nullify_grib(inputfile, outputfile, variables):
    """"
    Set to zero a variable in a GRIB file.
    This is done by unpacking the GRIB file, setting it to zero with CDO and then repacking it
    """ 

    print(f"Nullifying variable {variables} in GRIB file {inputfile}")
    
    if os.path.exists(inputfile):

        varlist=','.join(variables)
        singlefile = cdo.selname(varlist, input=inputfile, options="--eccodes")
        tempfile = cdo.mulc(0, input=singlefile, options="--eccodes")

        replace_field(inputfile, tempfile, outputfile, variables)
    else: 
        print(f'{inputfile} does not exist!')


def modify_grib(inputfile, outputfile, myfunction, spectral=False, **kwargs):
    """
    Modify a GRIB file using a function.
    Unpack grib1 and grib2, convert them to gaussian regular, 
    apply the function and the convert them back to grib1 and grib2.
    """

    # Unpack the GRIB file
    grib1, grib2 = unpack_grib_file(inputfile, "tmp")
    for file in [grib1, grib2]:
    
        if os.path.exists(file):
            print(f"Converting to netcdf file {file}")

            # Convert to netcdf: if spectral use sp2gpl, else use setgridtype
            if spectral:
                netcdf = cdo.sp2gpl(input=file, options=NC4)
            else:
                netcdf = cdo.setgridtype("regular", input=file, options=NC4)
            print(f"Modifying GRIB file {file} using function {myfunction.__name__}")

            # open the netcdf and modify it
            field = xr.open_dataset(netcdf,  engine="netcdf4")
            field = myfunction(field, **kwargs)
            
            # Save to a temporary file and remove
            with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
                temp_path = tmpfile.name
            field.to_netcdf(temp_path)
            shutil.move(temp_path, netcdf)
            
            print(f"Converting back to GRIB file {file}")
            grib = GRIB1 if file == grib1 else GRIB2
            if spectral:
                cdo.gp2spl(input=netcdf, output=file, options=grib)
            else:
                cdo.remapnn(inputfile, input=netcdf, output=file, options=grib)

    repack_grib_file(grib1, grib2, outputfile, clean=True)

def modify_single_grib(inputfile, outputfile, variables, myfunction, spectral=False, **kwargs):
    """
    Modify a GRIB file using a function.
    Unpack grib1 and grib2, convert them to gaussian regular, 
    apply the function and the convert them back to grib1 and grib2.
    """

    # Unpack the GRIB file
    
    if os.path.exists(inputfile):
        print(f"Converting to netcdf file {inputfile}")

        varlist=','.join(variables)
        singlefile = cdo.selname(varlist, input=inputfile, options="--eccodes")
        grib_version = detect_grib_version(singlefile)

        # Convert to netcdf: if spectral use sp2gpl, else use setgridtype
        if spectral:
            netcdf = cdo.sp2gpl(input=singlefile, options=NC4)
        else:
            netcdf = cdo.setgridtype("regular", input=singlefile, options=NC4)
        print(f"Modifying GRIB file {inputfile} using function {myfunction.__name__}")

        # open the netcdf and modify it
        field = xr.open_dataset(netcdf,  engine="netcdf4", decode_times=False)
        field = myfunction(field, var=variables, **kwargs)
        
        # Save to a temporary file and remove
        with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
            temp_path = tmpfile.name
        field.to_netcdf(temp_path)
        shutil.move(temp_path, netcdf)
        
        print(f"Converting back to GRIB file {singlefile}")
        grib = GRIB1 if grib_version==1 else GRIB2
        if spectral:
            cdo.gp2spl(input=netcdf, output=singlefile, options=grib)
        else:
            cdo.remapnn(inputfile, input=netcdf, output=singlefile, options=grib)

        replace_field(inputfile, singlefile, outputfile, variables)
    else: 
        print(f'{inputfile} does not exist!')

def truncate_grib_file(inputfile, outputfile, variables, orig=63, trunc=1):
    """
    Truncate the GRIB file to a specific size.
    """
    varlist=','.join(variables)
    print(varlist)
    trunc = cdo.sp2sp(str(trunc), input=f"-selname,{varlist} {inputfile}", options=NC4)
    compact = cdo.sp2sp(str(orig), input=trunc, options=GRIB2)
    where_expr = ",".join([f"shortName!={v}" for v in variables])
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
            temp_path = tmpfile.name
    subprocess.run(["grib_copy", "-w", where_expr, inputfile, temp_path], check=True)
    subprocess.run(["grib_copy", compact, temp_path, outputfile], check=True)


def detect_grib_version(filepath):
    """Detect the GRIB version of a file."""
    with open(filepath, "rb") as f:
        header = f.read(8)  # read enough bytes
        if header[:4] != b"GRIB":
            return None  # Not a GRIB file
        return header[7]  # Edition number (1 or 2)

def replace_field(inputfile, singlefile, outputfile, variable):
    """
    Replace a field in a GRIB file using grib_copy.
    """

    # allow for replacament
    if inputfile == outputfile:
        shutil.move(inputfile, "tmp.grib")
        inputfile = "tmp.grib"

    if os.path.exists(outputfile):
        os.remove(outputfile)
    if os.path.exists("filtered.grib"):
        os.remove("filtered.grib")
    if isinstance(variable, str):
        variable = [variable]
    where_expr = ",".join([f"shortName!={v}" for v in variable])
    subprocess.run([
        "grib_copy", "-w", where_expr,  # condition: where shortName is NOT t
        inputfile, "filtered.grib"], check=True)
    if os.path.exists("filtered.grib"):
        subprocess.run(["grib_copy", "filtered.grib", singlefile, outputfile], check=True)
        os.remove("filtered.grib")
    else:
        shutil.copyfile(singlefile, outputfile)
    os.remove(singlefile)
    

def unpack_grib_file(inputfile, tmpfile):
    """
    Unpack a GRIB file using grib_copy.
    This is used to split the GRIB messages into GRIB1 and GRIB2.
    """
    print(f"Unpacking GRIB file {inputfile} into {tmpfile}_grib1 and {tmpfile}_grib2")
    subprocess.run(["grib_copy", "-w", "edition=1", inputfile, f"{tmpfile}_grib1"], check=True)
    subprocess.run(["grib_copy", "-w", "edition=2", inputfile, f"{tmpfile}_grib2"], check=True)
    return f"{tmpfile}_grib1", f"{tmpfile}_grib2"

def repack_grib_file(grib1, grib2, outputfile, clean=True):

    """
    Repack a GRIB file using grib_copy.
    This is used to merge the GRIB1 and GRIB2 files back together.
    """
    if os.path.exists(outputfile):
        os.remove(outputfile)
    if not os.path.exists(grib1):
        print(f"Repacking GRIB file {grib2} into {outputfile}")
        subprocess.run(["grib_copy", grib2, outputfile], check=True)
    elif not os.path.exists(grib2):
        print(f"Repacking GRIB file {grib1} into {outputfile}")
        subprocess.run(["grib_copy", grib1, outputfile], check=True)
    else:
        print(f"Repacking GRIB file {grib1} and {grib2} into {outputfile}")
        subprocess.run(["grib_copy", grib1, grib2, outputfile], check=True)   

    # cleanup
    if clean:
        for file in [grib1, grib2]:
            if os.path.exists(file):
                os.remove(file)

def modify_value(field, var, newvalue):
    """
    Modify a field in the GRIB file setting it to a constant
    """
    for v in var:
        if v in field.variables:
            print(f"Modifying variable {v} in the field")
            field[v].data = np.full(field[v].shape, newvalue)
    return field

def replace_value(field, var, newfield):
    """
    Replace the field in a dataset with a dataarray
    """
    if 'time' not in newfield.dims:
        newfield = newfield.expand_dims('time', axis=0)

    newfield = newfield.transpose('time', 'lat', 'lon')

    for v in var:
        if v in field.variables:
            print(f"Replacing variable {v} in the field")
            field[v].data = newfield.data
    return field


def ecmwf_grid(kind):
    """Get the info on the grid to find the right ECMWF file"""
    
    ecmwf_name = {
        'L': 'l_2',
        'CO': '_4',
        'Q': '_2'
    }

    return ecmwf_name[kind.upper()]

def extract_grid_info(string):
    """Extract grid info from a string"""
    string = string.upper()
    pattern = r'T(CO|L)(\d+)L(\d+)'
    match = re.match(pattern, string)
    if match:
        grid_type = match.group(1)
        spectral = int(match.group(2))
        num_levels = int(match.group(3))
        return grid_type, spectral, num_levels
    
    return None
    
def spectral2gaussian(spectral, kind):
    """Convert spectral resolution to gaussian"""
    if kind.upper() == "CO":
        return int(spectral) + 1
    if kind == "L":
        return int((int(spectral) + 1) / 2)

    raise ValueError("Unknown grid type")


def albedo(field: xr.Dataset, var=None, newfield=None, **kwargs):
    if newfield is None:
        raise ValueError("You must provide `newfield` (land-sea mask).")

    print("Applying land–sea mask and zonal reconstruction...")

    # Broadcast land-sea mask to match field dims
    lsm_exp = newfield.broadcast_like(field)
    land_mask = lsm_exp > 0.5

    reconstructed = xr.Dataset()

    for v in field.data_vars:
        da = field[v]
        masked = da.where(land_mask)

        # Sort latitude ascending
        lat_ascending = da["lat"].values
        flip = False
        if lat_ascending[0] > lat_ascending[-1]:
            da = da.sortby("lat")
            masked = masked.sortby("lat")
            flip = True

        # Zonal mean
        zonal_mean = masked.mean(dim="lon", skipna=True)

        # Fill missing latitudes
        zonal_mean_filled = zonal_mean.interpolate_na(dim="lat", method="nearest")

        # Broadcast back to original shape
        da_recon = zonal_mean_filled.broadcast_like(da)

        # Apply polar band filling
        lat = da["lat"]
        band_s_mask = (lat < -52)
        if band_s_mask.any():
            mean_s = da_recon.sel(lat=lat.where((lat >= -52) & (lat <= -46), drop=True)).mean(dim=("lat","lon"), skipna=True)
            da_recon = da_recon.where(~band_s_mask, other=mean_s)

        band_n_mask = (lat > 75)
        if band_n_mask.any():
            mean_n = da_recon.sel(lat=lat.where((lat >= 70) & (lat <= 75), drop=True)).mean(dim=("lat","lon"), skipna=True)
            da_recon = da_recon.where(~band_n_mask, other=mean_n)

        # Flip back if needed
        if flip:
            da_recon = da_recon.sortby("lat", ascending=False)

        reconstructed[v] = da_recon

    print("Albedo modification complete, coordinates preserved for CDO.")
    return reconstructed

def eocene_mask(field, var, landsea=None, **kwargs):
    """
    Apply the Eocene land-sea mask to climate fields.
    Only modifies the mask grid (landsea), never the field grid.
    """

    if landsea is None:
        raise ValueError("Eocene land-sea mask must be provided as `landsea`.")

    # --- Extract data variable if needed ---
    if isinstance(landsea, xr.Dataset):
        mask_var = list(landsea.data_vars)[0]
        landsea = landsea[mask_var]

    # --- Ensure landsea has ascending latitude (xarray interp requires it) ---
    if not np.all(np.diff(landsea["lat"]) > 0):
        landsea = landsea.sortby("lat")

    # --- Interpolate the Eocene mask onto the ICMCL grid ---
    # Only the mask is modified here — field grid remains untouched.
    landsea_interp = landsea.interp(
        lat=field["lat"], lon=field["lon"], method="nearest"
    )

    # --- Create boolean mask ---
    mask_exp = (landsea_interp > 0.5).broadcast_like(field)

    result = field.copy()

    # --- Apply rules ---
    albedo_vars = ["al", "aluvp", "aluvd", "alnip", "alnid"]
    lai_vars = ["lai_lv", "lai_hv"]

    # Albedo rule: ocean = 0.05
    for v in albedo_vars:
        if v in result:
            result[v] = result[v].where(mask_exp, 0.05)

    # LAI rule: ocean = NaN
    for v in lai_vars:
        if v in result:
            result[v] = result[v].where(mask_exp)

    return result




