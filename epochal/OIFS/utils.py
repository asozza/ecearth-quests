"""Some utilities for OIFS grid definition"""
import re
import os
import tempfile
import shutil
import eccodes
import numpy as np
import xarray as xr
import xesmf as xe
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

        grib_version = detect_grib_version(inputfile)
        grib_opt = GRIB1 if grib_version == 1 else GRIB2


        varlist=','.join(variables)
        singlefile = cdo.selname(varlist, input=inputfile, options="--eccodes")
        
        # Convert to netcdf: if spectral use sp2gpl, else use setgridtype
        if spectral:
            netcdf = cdo.sp2gpl(input=singlefile, options=NC4)
        else:
            netcdf = cdo.setgridtype("regular", input=singlefile, options=NC4)
    #if spectral:
    #    netcdf = cdo.sp2gpl(input=inputfile, options=NC4)
    #else:
    #    netcdf = cdo.setgridtype("regular", input=inputfile, options=NC4)
    print(f"Modifying GRIB file {inputfile} using function {myfunction.__name__}")

        # open the netcdf and modify it
    field = xr.open_dataset(netcdf,  engine="netcdf4", decode_times=False)
    print(f"Applying {myfunction.__name__} to variables {variables}")
    field = myfunction(field, var=variables, **kwargs)
        
        # Save to a temporary file and remove
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        temp_path = tmpfile.name
    field.to_netcdf(temp_path)
    shutil.move(temp_path, netcdf)
        
    #with tempfile.NamedTemporaryFile(delete=False) as tmp:
    #        tmp_nc = tmp.name
    #field.to_netcdf(tmp_nc)

    print(f"Converting back to GRIB file {singlefile}")
        #grib = GRIB1 if grib_version==1 else GRIB2
    if spectral:
        cdo.gp2spl(input=netcdf, output=singlefile, options=grib_opt)
    #else:
    #    cdo.remapnn(inputfile, input=netcdf, output=singlefile, options=grib_opt)

    #replace_field(inputfile, singlefile, outputfile, variables)
    #else: 
        #print(f'{inputfile} does not exist!')
    #grib_opt = GRIB1 if grib_version == 1 else GRIB2
    #if spectral:
    #    cdo.gp2spl(input=tmp_nc, output=outputfile, options=grib_opt)
    else:
        cdo.copy(input=netcdf, output=singlefile, options="--eccodes")

    #print(f"Updated GRIB written to {outputfile}")
    replace_field(inputfile, singlefile, outputfile, variables)

def modify_single_grib_(
    inputfile,
    outputfile,
    variables,
    myfunction,
    spectral=False,
    **kwargs
):
    """
    Modify one or more variables in a GRIB file while preserving GRIB identity.
    Includes verbose debug output for tracing execution.
    """

    print("\n" + "=" * 80)
    print("START modify_single_grib_")
    print(f"Input GRIB   : {inputfile}")
    print(f"Output GRIB  : {outputfile}")
    print(f"Variables    : {variables}")
    print(f"Spectral     : {spectral}")
    print("=" * 80)

    if isinstance(variables, str):
        variables = [variables]

    if not os.path.exists(inputfile):
        raise FileNotFoundError(inputfile)

    # --------------------------------------------------
    # 0. Extract GRIB identity
    # --------------------------------------------------
    print("\n[0] Extracting GRIB metadata (paramId + shortName)")

    metadata = {}

    with open(inputfile, "rb") as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break

            try:
                shortName = eccodes.codes_get(gid, "shortName")
            except eccodes.KeyValueNotFoundError:
                shortName = None

            if shortName in variables:
                paramId = eccodes.codes_get(gid, "paramId")
                metadata[shortName] = {
                    "shortName": shortName,
                    "paramId": paramId,
                }
                print(f"  ✔ Found {shortName}: paramId={paramId}")

            eccodes.codes_release(gid)

    if not metadata:
        raise RuntimeError("None of the requested variables were found in the GRIB file")

    # --------------------------------------------------
    # 1. Extract requested fields
    # --------------------------------------------------
    varlist = ",".join(variables)
    print(f"\n[1] Extracting variable(s): {varlist}")
    single_grib = cdo.selname(varlist, input=inputfile, options="--eccodes")
    print(f"    → Temporary GRIB: {single_grib}")

    # --------------------------------------------------
    # 2. Convert GRIB → NetCDF
    # --------------------------------------------------
    print("\n[2] Converting GRIB → NetCDF")

    if spectral:
        print("    Using spectral conversion (sp2gpl)")
        nc_in = cdo.sp2gpl(input=single_grib, options=NC4)
    else:
        print("    Using grid-point conversion (setgridtype,regular)")
        nc_in = cdo.setgridtype("regular", input=single_grib, options=NC4)

    print(f"    → NetCDF file: {nc_in}")

    # --------------------------------------------------
    # 3. Modify using xarray
    # --------------------------------------------------
    print("\n[3] Opening NetCDF and applying modification function")

    ds = xr.open_dataset(nc_in, decode_times=False)
    print(f"    Variables in NetCDF: {list(ds.data_vars)}")

    ds = myfunction(ds, var=variables, **kwargs)
    print(f"    ✔ Modification function '{myfunction.__name__}' applied")

    # --------------------------------------------------
    # 3b. Re-inject GRIB metadata
    # --------------------------------------------------
    print("\n[3b] Re-injecting GRIB metadata into NetCDF variables")

    for v in variables:
        if v not in ds:
            print(f"    ⚠ Variable {v} not found in dataset after modification")
            continue

        meta = metadata[v]
        ds[v].attrs["GRIB_paramId"] = meta["paramId"]
        ds[v].attrs["GRIB_shortName"] = meta["shortName"]

        print(
            f"    ✔ {v}: GRIB_paramId={meta['paramId']}, "
            f"GRIB_shortName={meta['shortName']}"
        )

    # --------------------------------------------------
    # 4. Write modified NetCDF
    # --------------------------------------------------
    print("\n[4] Writing modified NetCDF to temporary file")

    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
        nc_out = tmp.name

    ds.to_netcdf(nc_out)
    ds.close()

    print(f"    → NetCDF written: {nc_out}")

    # --------------------------------------------------
    # 5. Convert NetCDF → GRIB
    # --------------------------------------------------
    print("\n[5] Converting NetCDF → GRIB")

    grib_version = detect_grib_version(inputfile)
    grib_opt = GRIB1 if grib_version == 1 else GRIB2

    print(f"    Detected GRIB edition: {grib_version}")
    print(f"    CDO options used     : {grib_opt}")

    if spectral:
        cdo.gp2spl(input=nc_out, output=single_grib, options=grib_opt)
    else:
        cdo.copy(input=nc_out, output=single_grib, options=grib_opt)

    print(f"    ✔ GRIB updated: {single_grib}")

    # --------------------------------------------------
    # 6. Merge modified field into full GRIB
    # --------------------------------------------------
    print("\n[6] Merging modified field back into original GRIB")

    replace_field(
        inputfile=inputfile,
        singlefile=single_grib,
        outputfile=outputfile,
        variable=variables,
    )

    print(f"    ✔ Field replaced in output GRIB")

    print("\n[FINAL] Forcing variable name with cdo setname")

    tmp = tempfile.NamedTemporaryFile(delete=False).name
    cdo.setname(",".join(variables), input=outputfile, output=tmp)
    shutil.move(tmp, outputfile)

    print("✔ Variable name forced successfully")

    # --------------------------------------------------
    # 7. Cleanup
    # --------------------------------------------------
    print("\n[7] Cleaning up temporary files")

    os.remove(nc_out)
    print(f"    Removed temporary NetCDF: {nc_out}")

    print("\n" + "=" * 80)
    print("END modify_single_grib_")
    print("=" * 80 + "\n")


def _ensure_paramid(ds, variables):
    """
    Attach ECMWF paramId metadata so CDO does not write codeXXX variables.
    Extend PARAMIDS as needed.
    """

    PARAMIDS = {
        "sdor": 160,   # standard deviation of sub-grid orography
        "lsm": 172,    # land-sea mask
        # add more ECMWF local params here if needed
    }

    for v in variables:
        if v in ds and v in PARAMIDS:
            ds[v].attrs["paramId"] = PARAMIDS[v]
            ds[v].attrs["shortName"] = v

    return ds



def modify_grib_file(inputfile, outputfile, variables, myfunction, spectral=False, **kwargs):

    grib_version = detect_grib_version(inputfile)
    grib_opt = GRIB1 if grib_version == 1 else GRIB2

    # 1️⃣ Convert FULL GRIB → NetCDF WITHOUT changing grid
    netcdf = cdo.copy(input=inputfile, options=NC4)

    # 2️⃣ Modify with xarray
    ds = xr.open_dataset(netcdf, decode_times=False)
    ds = myfunction(ds, var=variables, **kwargs)

    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        ds.to_netcdf(tmp.name)

        # 3️⃣ Convert back → GRIB, preserving structure
        if spectral:
            cdo.gp2spl(input=tmp.name, output=outputfile, options=grib_opt)
        else:
            cdo.copy(input=tmp.name, output=outputfile, options=grib_opt)

    print(f"Updated GRIB written to {outputfile}")

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
    Replace variables in a dataset with values from a DataArray or single-variable Dataset.
    """

    # Unwrap Dataset → DataArray if needed
    if isinstance(newfield, xr.Dataset):
        if len(newfield.data_vars) != 1:
            raise ValueError("newfield Dataset must contain exactly one variable")
        newdata = next(iter(newfield.data_vars.values()))
    else:
        newdata = newfield

    # Ensure time dimension exists
    if 'time' not in newdata.dims:
        newdata = newdata.expand_dims('time', axis=0)

    # Ensure correct dimension order
    newdata = newdata.transpose('time', 'lat', 'lon')

    for v in var:
        if v in field.variables:
            print(f"Replacing variable {v} in the field")

            # Safety check
            if field[v].shape != newdata.shape:
                raise ValueError(
                    f"Shape mismatch for {v}: "
                    f"{field[v].shape} vs {newdata.shape}"
                )

            field[v].data = newdata.data

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


def regrid_dataset(data, regrid_to_reference):
    """
    Regrid a DataArray to match the grid of another DataArray.
    """
    regridder = xe.Regridder(
        data, 
        regrid_to_reference, 
        method='bilinear'
    )
    return regridder(data)


def remap_to_field_grid(source_da, target_field, method="nn"):
    """
    Remap source_da onto the grid of target_field using CDO.
    """

    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as src, \
         tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tgt, \
         tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as grd, \
         tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as out:

        # Write files
        source_da.to_netcdf(src.name)
        target_field.to_netcdf(tgt.name)

        # Extract GRID DESCRIPTION (this is the key step)
        cdo.griddes(input=tgt.name, output=grd.name)

        # Remap
        if method == "nn":
            cdo.remapnn(grd.name, input=src.name, output=out.name)
        elif method == "bil":
            cdo.remapbil(grd.name, input=src.name, output=out.name)
        else:
            raise ValueError("Unknown remap method")

        ds = xr.open_dataset(out.name)
        return list(ds.data_vars.values())[0]


def get_grib_metadata(inputfile, variables):
    """
    Return a dict mapping variable name → GRIB metadata (shortName, paramId, discipline, category, number)
    """
    metadata = {}
    for v in variables:
        with open(inputfile, "rb") as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    name = eccodes.codes_get(gid, "shortName")
                except eccodes.KeyValueNotFoundError:
                    name = None

                if name == v:
                    meta_dict = {}
                    for key in ["shortName", "paramId", "discipline", "category", "number"]:
                        try:
                            meta_dict[key] = eccodes.codes_get(gid, key)
                        except eccodes.KeyValueNotFoundError:
                            meta_dict[key] = None
                    metadata[v] = meta_dict

                eccodes.codes_release(gid)
    return metadata


