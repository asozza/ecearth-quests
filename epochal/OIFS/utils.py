"""Some utilities for OIFS grid definition"""
import re
import os
import tempfile
import shutil
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
        cdo.copy(input=tempfile, output="nulify.nc")
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
        cdo.copy(input=singlefile, output='singlefile.nc')

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
        #os.remove("filtered.grib")
    else:
        shutil.copyfile(singlefile, outputfile)
    #os.remove(singlefile)
    

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