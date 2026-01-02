#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Command line tool to modify the NEMO restart files from a EC-Eart4 experiment, 
given a specific experiment name and time leg. 

Needed modules:
# module load intel/2021.4.0 intel-mkl/19.0.5 prgenv/intel hdf5 netcdf4
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/apps/netcdf4/4.9.1/INTEL/2021.4/lib:/usr/local/apps/hdf5/1.12.2/INTEL/2021.4/lib

Authors: Alessandro Sozza (CNR-ISAC)
Date: Nov 2023
"""

import os
import glob
import time
import yaml
import cftime
import shutil
import logging
import argparse
import platform
import subprocess
import datetime

import numpy as np
from cdo import Cdo
import xarray as xr
import netCDF4 as nc

# Initialize CDO
cdo = Cdo()

# routines / sections:
# [1] time utils, [2] file utils
# [3] field reader, [4] restart rebuild & reader
# [5] field averaging and manipulation
# [6] cdo, [7] EOF, [8] forecaster
# [9] parser

# Define global varlists for each variable
varlists = {
    'thetao': ['tn', 'tb'],
    'so': ['sn', 'sb'],
    'zos': ['sshn', 'sshb'],
    'uo': ['un', 'ub'],
    'vo': ['vn', 'vb']
}

# define minimal catalogue
catalogue = {
    'thetao' : {'dim': '3D', 'grid': 'T', 'units': 'degC', 'long_name': 'Temperature'},
    'so': {'dim': '3D', 'grid': 'T', 'units': 'PSU', 'long_name': 'Salinity'}
}

# define global paths
base_path = "/ec/res4/scratch/itas/ece4"
src_path = "/ec/res4/hpcperm/itas/src/gitlab/ecearth4-v4.1.3/sources/nemo-4.2/tools/REBUILD_NEMO"
data_path = "/lus/h2resw01/hpcperm/ccpd/ECE4-DATA"

def folders(expname):
    """ List of global paths dependent on expname """

    dirs = {
        'exp': os.path.join(base_path, expname),
        'nemo': os.path.join(base_path, expname, "output", "nemo"),
        'oifs': os.path.join(base_path, expname, "output", "oifs"),
        'restart': os.path.join(base_path, expname, "restart"),
        'log': os.path.join(base_path, expname, "log"),
        'tmp': os.path.join(base_path, expname, "tmp"),
        'post': os.path.join(base_path, expname, "post"),
        'domain': os.path.join(data_path, "nemo", "domain")
    }

    # Create 'post' & 'tmp' folder if it doesn't exist
    if not os.path.exists(dirs['post']):
        os.makedirs(dirs['post'])

    if not os.path.exists(dirs['tmp']):
        os.makedirs(dirs['tmp'])

    return dirs

##########################################################################################
# [1] time utils

def get_epoch(date):
    """ Get epoch from date """

    return time.mktime(date.timetuple())

def get_year_fraction(date):
    """ Transform date into year fraction """

    start_of_year = datetime.datetime(date.year,1,1,0,0,0)
    end_of_year = datetime.datetime(date.year+1,1,1,0,0,0)
    year_elapsed = get_epoch(date) - get_epoch(start_of_year)
    year_duration = get_epoch(end_of_year) - get_epoch(start_of_year)
    Frac = year_elapsed/year_duration

    return  date.year + Frac

def get_decimal_year(date):
    """ Get decimal year from year fraction """

    return [get_year_fraction(d) for d in date]

def _forecast_xarray(year):
    """Get the xarray for the forecast time"""

    fdate = cftime.DatetimeGregorian(year, 1, 1, 0, 0, 0, has_year_zero=False)
    
    xf = xr.DataArray(data = np.array([fdate]), dims = ['time'], coords = {'time': np.array([fdate])},
                      attrs = {'stardand_name': 'time', 'long_name': 'Time axis', 'bounds': 'time_counter_bnds', 'axis': 'T'})

    return xf

##########################################################################################
# [2] file utils

def delete_attrs(file):
    """ Delete attributes """

    # Open the dataset in 'r+' mode to allow modifications
    with nc.Dataset(file, 'r+') as dataset:
        # Iterate over all variables in the dataset
        for var_name, variable in dataset.variables.items():
            # Get all attribute names for the variable
            attr_names = list(variable.ncattrs())
            # Delete each attribute
            for attr in attr_names:
                variable.delncattr(attr)  # Correct method to delete an attribute

    return None

def remove_existing_file(filename):
    """ Remove a file if it exists """

    try:
        os.remove(filename)
        logging.info(f"File {filename} successfully removed.")
    except FileNotFoundError:
        logging.info(f"File {filename} not found.")

def remove_existing_filelist(filename):
    """ Remove all files matching the filename pattern """

    pattern = os.path.join(filename + '*.nc')
    files = glob.glob(pattern)
    try:
        for file in files:
            os.remove(file)
            logging.info(f"File {file} successfully removed.")
    except FileNotFoundError:
        logging.error(f"File {file} not found.")

def get_nemo_timestep(filename):
    """ Get timestep from a NEMO restart file """

    return os.path.basename(filename).split('_')[1]

def cleanup_eofs(expname, varname, endyear, window, yearleap):
    """Clean up routine for eofs"""

    dirs = folders(expname)
    dirs['tmp'] = os.path.join(dirs['tmp'], f'Y{endyear}-W{window}-L{yearleap}')

    filename = os.path.join(dirs['tmp'], f"{varname}.nc")
    remove_existing_file(filename)
    filename = os.path.join(dirs['tmp'], f"{varname}_anomaly.nc")
    remove_existing_file(filename)
    filename = os.path.join(dirs['tmp'], f"{varname}_pattern.nc")
    remove_existing_file(filename)
    filename = os.path.join(dirs['tmp'], f"{varname}_variance.nc")
    remove_existing_file(filename)
    timeseries = os.path.join(dirs['tmp'], f"{varname}_series_")
    remove_existing_filelist(timeseries)

    return None

##########################################################################################
# [3] reader routines

def _nemodict(grid, freq):
    """ 
    Nemodict: Dictionary of NEMO output fields
    
    Args: 
    grid: grid name [T, U, V, W]
    freq: output frequency [1m, 1y, ...]

    """

    gridlist = ["T", "U", "V"]    
    if grid in gridlist:
        grid_lower = grid.lower()
        return {
            grid: {
                "preproc": preproc_nemo,
                "format": f"oce_{freq}_{grid}",
                "x_grid": f"x_grid_{grid}",
                "y_grid": f"y_grid_{grid}",
                "nav_lat": f"nav_lat_grid_{grid}",
                "nav_lon": f"nav_lon_grid_{grid}",
                "depth": f"depth{grid_lower}",
                "x_grid_inner": f"x_grid_{grid}_inner",
                "y_grid_inner": f"y_grid_{grid}_inner"
            }
        }
    elif grid == "W":
        grid_lower = grid.lower()
        return {
            "W": {
                "preproc": preproc_nemo,
                "format": f"oce_{freq}_{grid}",
                "nav_lat": "nav_lat",
                "nav_lon": "nav_lon",
                "depth": f"depth{grid_lower}"
            }
        }
    elif grid == "ice":
        return {
            "ice": {
                "preproc": preproc_nemo_ice,
                "format": f"ice_{freq}"
            }
        }
    else:
        raise ValueError(f"Unsupported grid type: {grid}")

def preproc_nemo(data, grid):
    """ 
    General preprocessing routine for NEMO data based on grid type
    
    Args: 
    data: dataset
    grid: gridname [T, U, V, W]

    """
    
    grid_mappings = _nemodict(grid, None)[grid]  # None for freq as it is not used here

    if grid != 'W':
        data = data.rename_dims({grid_mappings["x_grid"]: 'x', grid_mappings["y_grid"]: 'y'})
        data = data.swap_dims({grid_mappings["x_grid_inner"]: 'x', grid_mappings["y_grid_inner"]: 'y'})

    data = data.rename({
        grid_mappings["nav_lat"]: 'lat', 
        grid_mappings["nav_lon"]: 'lon', 
        grid_mappings["depth"]: 'z', 
        'time_counter': 'time'
    })
    data = data.drop_vars(['time_centered'], errors='ignore')
    data = data.drop_dims(['axis_nbounds'], errors='ignore')

    return data

def preproc_nemo_ice(data):
    """Preprocessing routine for NEMO for ice"""

    data = data.rename({'time_counter': 'time'})
    
    return data

def reader_nemo(expname, startyear, endyear, grid="T", freq="1m"):
    """ 
    reader_nemo: function to read NEMO data 
    
    Args:
    expname: experiment name
    startyear,endyear: time window
    grid: grid name [T, U, V, W]
    frequency: output frequency [1m, 1y, ...]

    """

    dirs = folders(expname)
    dict = _nemodict(grid, freq)

    filelist = []
    available_years = []
    for year in range(startyear, endyear + 1):
        pattern = os.path.join(dirs['nemo'], f"{expname}_{dict[grid]['format']}_{year}-{year}.nc")
        matching_files = glob.glob(pattern)
        if matching_files:
            filelist.extend(matching_files)
            available_years.append(year)

    if not filelist:
        raise FileNotFoundError(f"No data files found for the specified range {startyear}-{endyear}.")

    # Log a warning if some years are missing
    if available_years:
        actual_startyear = min(available_years)
        actual_endyear = max(available_years)
        if actual_startyear > startyear or actual_endyear < endyear:
            logging.warning(f"Data available only in the range {actual_startyear}-{actual_endyear}.")
        else:
            logging.info(f"Data available in the range {startyear}-{endyear}.")
    else:
        raise FileNotFoundError("No data files found within the specified range.")

    #logging.info('Files to be loaded %s', filelist)
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    data = xr.open_mfdataset(filelist, preprocess=lambda d: dict[grid]["preproc"](d, grid), decode_times=time_coder, data_vars="all")

    return data


def reader_nemo_field(expname, startyear, endyear, varname, freq="1m"):
    """ 
    reader_nemo_field: function to read NEMO field 
    
    Args:
    expname: experiment name
    startyear,endyear: time window
    varname: variable name

    """

    info = catalogue[varname]

    # Check for 'dependencies' if dealing with derived variable 
    if 'dependencies' in info: 
        field = {}
        for grid, var in zip(info['grid'], info['dependencies']):
            data = reader_nemo(expname=expname, startyear=startyear, endyear=endyear, grid=grid)
            field[var] = data[var]
            if 'preprocessing' in info and var in info['preprocessing']:
                field[var] = info['preprocessing'][var](field[var])
        data = info['operation'](*[field[var] for var in info['dependencies']])
    else:
        data = reader_nemo(expname=expname, startyear=startyear, endyear=endyear, grid=info['grid'], freq=freq)
        data = data[[varname]]

    return data

##########################################################################################
# [4] restart routines

def rebuilder(expname, leg):
    """Function to rebuild NEMO restart files for a given experiment and time leg.

    This function processes the restart files for the specified experiment and time leg,
    creates necessary symbolic links, and runs the rebuild_nemo executable to rebuild
    the NEMO restart files. It handles both 'restart' and 'restart_ice' file types.

    Args:
        expname (str): The name of the experiment.
        leg (str): The time leg of the experiment.

    Raises:
        subprocess.CalledProcessError: If the rebuild_nemo command fails.
    """

    dirs = folders(expname)
    
    os.makedirs(os.path.join(dirs['tmp'], str(leg).zfill(3)), exist_ok=True)

    rebuild_exe = os.path.join(src_path, "rebuild_nemo")
  
    for kind in ['restart', 'restart_ice']:
        print(' Processing ' + kind)
        flist = glob.glob(os.path.join(dirs['restart'], str(leg).zfill(3), expname + '*_' + kind + '_????.nc'))
        tstep = get_nemo_timestep(flist[0])

        for filename in flist:
            destination_path = os.path.join(dirs['tmp'], str(leg).zfill(3), os.path.basename(filename))
            try:
                os.symlink(filename, destination_path)
            except FileExistsError:
                pass

        rebuild_command = [rebuild_exe, "-m", os.path.join(dirs['tmp'], str(leg).zfill(3), expname + "_" + tstep + "_" + kind ), str(len(flist))]
        try:
            subprocess.run(rebuild_command, stderr=subprocess.PIPE, text=True, check=True)
            for file in glob.glob('nam_rebuld_*'):
                os.remove(file)
        except subprocess.CalledProcessError as e:
            error_message = e.stderr
            print(error_message)

        for filename in flist:
            destination_path = os.path.join(dirs['tmp'], str(leg).zfill(3), os.path.basename(filename))
            os.remove(destination_path)

    # copy restart
    tstep = get_nemo_timestep(glob.glob(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '*_restart.nc'))[0])
    shutil.copy(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + tstep + '_restart.nc'), os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart.nc'))
    shutil.copy(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + tstep + '_restart_ice.nc'), os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart_ice.nc'))

    # delete temporary files
    flist = glob.glob('nam_rebuild*')
    for file in flist:
        os.remove(file)

    return None

def reader_rebuilt(expname, startleg, endleg):
    """ Read rebuilt NEMO restart files """

    dirs = folders(expname)
    
    filelist = []
    for leg in range(startleg,endleg+1):
        pattern = os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '*_restart.nc')
        #pattern = os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart.nc')        
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    logging.info(' File to be loaded %s', filelist)
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    data = xr.open_mfdataset(filelist, engine="netcdf4", decode_times=time_coder, data_vars="all")

    return data


##########################################################################################
# [5] field averaging & manipulation routines

def window_mean(ds0, window):
    """ Compute the mean of a window of months in a dataset """

    time = ds0.time
    years = time.dt.year
    months = time.dt.month

    # Determine if the window is contained in one year or crosses two years
    min_month, max_month = window[0], window[-1]
    cross_year = min_month > max_month  # True if the window crosses two years (e.g. NDJF)

    grouped_data = {}
    new_time = []

    for year in set(years.values):
        if cross_year:
            # If the window crosses two years, months are divided in (year,year+1)
            mask = ((years == year) & (months >= min_month)) | \
                    ((years == year + 1) & (months <= max_month))
            ref_year = year  # Starting year of the window
        else:
            # If window is entirely in the same year
            mask = (years == year) & (months.isin(window))
            ref_year = year

        selected_ds = ds0.sel(time=mask)

        if selected_ds.time.size > 0:
            # Compute time mean for the window
            mean_value = selected_ds.mean(dim='time', keep_attrs=True)
            
            # Compute time centroid of the window
            time_values = selected_ds.time.values
            centroid_num = sum(cftime.date2num(t, units="seconds since 1990-01-01", calendar=t.calendar) for t in time_values) / len(time_values)
            centroid = cftime.num2date(centroid_num, units="seconds since 1990-01-01", calendar=time_values[0].calendar)

            # Use centroid as key (in cftime not string format)
            grouped_data[centroid] = mean_value
            new_time.append(centroid)  # maintain cftime

    # Create new dataset with annual means and new temporal coordinates
    new_ds = xr.concat(list(grouped_data.values()), dim='time')
    new_ds = new_ds.assign_coords(time=new_time)  # Assign new dates in cftime
    new_ds = new_ds.sortby('time')

    return new_ds

def constraints_for_fields(data):
    """ 
    Check and apply constraints to variables 
    
    U < 10 m/s, |ssh| < 20 m, S in [0,100] psu, T > -2.5 degC
    """ 

    # for horizontal velocity (u,v)
    for var in ['uo', 'vo']:
        if var in data:
            data[var] = data[var].clip(-10.0,10.0) # Ensure U < 10 m/s
    
    # for sea surface height (ssh)
    for var in ['zos']:
        if var in data:
            data[var] = data[var].clip(-20,20)
    
    # for salinity
    for var in ['so']:
        if var in data:
            data[var] = data[var].clip(0, 100)  # Ensure S in [0, 100]
    
    # for temperature
    for var in ['thetao']:
        if var in data:
            data[var] = data[var].clip(-2.5, 32.0) # Ensure T in [-2.5, 32]

    return data

##########################################################################################
# [6] cdo routines

def merge(expname, varname, startyear, endyear, window, yearleap):
    """ CDO command to merge data """
    
    dirs = folders(expname)

    dirs['tmp'] = os.path.join(dirs['tmp'], f'Y{endyear}-W{window}-L{yearleap}')
    os.makedirs(dirs['tmp'], exist_ok=True)
        
    logging.info(f"pymerge start/end year: {startyear}-{endyear}") 
    data = reader_nemo_field(expname=expname, startyear=startyear, endyear=endyear, varname=varname)

    data_winter = window_mean(data, [12, 1])
    data_winter = data_winter.sel(time=data_winter['time.month'].isin(1))
    data_winter.to_netcdf(os.path.join(dirs['tmp'], f"{varname}.nc"))
    
    data_winter = data_winter - data_winter.mean(dim='time')
    data_winter.to_netcdf(os.path.join(dirs['tmp'], f"{varname}_anomaly.nc"))

    return None


def get_eofs(expname, varname, endyear, window, neofs, yearleap):
    """Compute EOF using the CDO Python package with error handling."""

    info = catalogue[varname]

    # Get the directories and file paths
    dirs = folders(expname)
    dirs['tmp'] = os.path.join(dirs['tmp'], f'Y{endyear}-W{window}-L{yearleap}')

    logging.info(f"Endyear: {endyear}")
    logging.info(f"Forecast window: {window}")
    logging.info(f"Number of EOFs: {neofs}")

    # Define file paths for anomaly, covariance, and pattern output files
    flda = os.path.join(dirs['tmp'], f"{varname}_anomaly.nc")
    fldcov = os.path.join(dirs['tmp'], f"{varname}_variance.nc")
    fldpat = os.path.join(dirs['tmp'], f"{varname}_pattern.nc")
    
    # Remove existing output files if they already exist
    remove_existing_file(fldcov)
    remove_existing_file(fldpat)

    logging.info(f"Computing EOFs for variable {varname} with window size {neofs}.")
    
    # CDO command to compute EOFs covariance and pattern
    if info['dim'] == '3D':
        logging.info(f"Compute 3D EOFs")
        cdo.run(f"eof3d,{neofs} {flda} {fldcov} {fldpat}")
    elif info['dim'] == '2D':   
        logging.info(f"Compute 2D EOFs")        
        cdo.run(f"eof,{neofs} {flda} {fldcov} {fldpat}")

    # Define timeseries output file pattern
    timeseries = os.path.join(dirs['tmp'], f"{varname}_series_")
    remove_existing_filelist(timeseries)

    # Compute EOF coefficients (timeseries)
    logging.info(f"Computing EOF coefficients for {varname}.")
    if info['dim'] == '3D':
        cdo.run(f"eofcoeff3d {fldpat} {flda} {timeseries}")
    elif info['dim'] == '2D':
        cdo.run(f"eofcoeff {fldpat} {flda} {timeseries}")

    logging.info(f"EOF computation completed successfully: {fldcov}, {fldpat}, {timeseries}")
    
    return None

##########################################################################################
# [7] EOF xarray routines

def _grid_mapping(grid):
    """
    Generate grid-specific renaming conventions for NEMO data.

    Args:
    - grid: str
        Grid type. Choose from 'T', 'U', 'V', 'W'.

    Returns:
    - Dictionary with renaming conventions for the specified grid.
    """

    grid_lower = grid.lower()
    if grid in {'T', 'U', 'V', 'W'}:
        return {
            "x": f"x_grid_{grid}",
            "y": f"y_grid_{grid}",
            "lat": f"nav_lat_grid_{grid}",
            "lon": f"nav_lon_grid_{grid}",
            "z": f"depth{grid_lower}"
        }
    else:
        raise ValueError(f"Unsupported grid type: {grid}")


def process_data(data, ftype, dim='3D', grid='T'):
    """
    Function for preprocessing or postprocessing 2D/3D data.

    Parameters:
    - data: xarray.DataArray or xarray.Dataset
        Input data to be processed.
    - ftype: str
        Operation type. Choose from:
        'pattern' - Preprocessing for EOF pattern.
        'series' - Preprocessing for EOF timeseries.
        'post' - Post-processing for variable field.
    - dim: str, optional (default='3D')
        Dimensionality of the dataset. Choose from '2D' or '3D'.
    - grid: str, optional (default='T')
        Grid type. Choose from 'T', 'U', 'V', 'W'.
        
    Returns:
    - Processed data (xarray.DataArray or xarray.Dataset)
    """

    if dim not in {'2D', '3D'}:
        raise ValueError(f"Invalid dim '{dim}'. Choose '2D' or '3D'.")

    if grid not in {'T', 'U', 'V', 'W'}:
        raise ValueError(f"Invalid grid '{grid}'. Choose from 'T', 'U', 'V', 'W'.")

    grid_mappings = _grid_mapping(grid)

    if ftype == 'pattern':
        # Preprocessing routine for EOF pattern
        data = data.rename_dims({grid_mappings['x']: 'x', grid_mappings['y']: 'y'})
        data = data.rename({grid_mappings['lat']: 'lat', grid_mappings['lon']: 'lon'})
        data = data.rename({'time_counter': 'time'})
        data = data.drop_vars({'time_counter_bnds'}, errors='ignore')
        if dim == '3D':
            data = data.rename({grid_mappings['z']: 'z'})
            data = data.drop_vars({f"{grid_mappings['z']}_bnds"}, errors='ignore')


    elif ftype == 'series':
        # Preprocessing routine for EOF timeseries
        #data = data.rename({'time_counter': 'time'})
        data = data.isel(lon=0, lat=0)
        data = data.drop_vars({'time_counter_bnds', 'lon', 'lat'}, errors='ignore')
        if dim == '3D':
            data = data.isel(zaxis_Reduced=0)
            data = data.drop_vars({'zaxis_Reduced'}, errors='ignore')

    elif ftype == 'post':
        # Post-processing routine for variable field
        data = data.rename({'time': 'time_counter'})
        data = data.rename_dims({'x': grid_mappings['x'], 'y': grid_mappings['y']})
        data = data.rename({'lat': grid_mappings['lat'], 'lon': grid_mappings['lon']})        
        if dim == '3D':
            data = data.rename({'z': grid_mappings['z']})

    else:
        raise ValueError(f"Invalid mode '{ftype}'. Choose from 'pattern', 'series', or 'post'.")

    return data


def project_eofs(expname, varname, endyear, window, neofs, yearleap, mode='full'):
    """ 
    Function to project a field in the future using EOFs and linear regression. 
    
    Different options are available:
    - projection using the full set of EOFs (mode=full)
    - projection using only the first EOF (mode=first)
    - projection using EOFs up to a percentage of the full (mode=frac)
    - reconstruction of the original field using EOFs (mode=reco)
    
    """

    endleg = endyear - 1990 + 2
    startleg = endleg - window + 1
    startyear = 1990 + startleg - 2
    targetyear = endyear + yearleap

    logging.info(f"Start/end year: {startyear}-{endyear}")
    logging.info(f"Regression window: {window}")

    dirs = folders(expname)
    dirs['tmp'] = os.path.join(dirs['tmp'], f'Y{endyear}-W{window}-L{yearleap}')

    info = catalogue[varname]

    # forecast target year
    xf = _forecast_xarray(targetyear)

    # prepare patterns for EOFs
    filename = os.path.join(dirs['tmp'], f"{varname}_pattern.nc")
    print(filename)
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    pattern = xr.open_dataset(filename, decode_times=time_coder)
    field = pattern.isel(time=0)*0


    # Full set of EOFs
    if mode == 'full':

        logging.info(f"Number of EOFs: {neofs}")
        for i in range(neofs):
            filename = os.path.join(dirs['tmp'], f"{varname}_series_{str(i).zfill(5)}.nc")
            logging.info(f"Reading filename: {filename}")
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            timeseries = xr.open_mfdataset(filename, decode_times=time_coder, preprocess=lambda data: process_data(data, ftype='series', dim=info['dim'], grid=info['grid']))
            coeffs = timeseries.polyfit(dim='time', deg=1, skipna=True)
            theta = xr.polyval(xf, coeffs[f"{varname}_polyfit_coefficients"])
            basis = pattern.isel(time=i)
            field += theta * basis
        
        field['time'] = xf

        # add mean (TO BE CHANGED!!!!!!!!!)
        logging.info(f"Adding mean trend of {varname}.nc")
        filename = os.path.join(dirs['tmp'], f"{varname}.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        data = xr.open_dataset(filename, decode_times=time_coder)
        #data = xr.open_mfdataset(filename) # use_cftime=True, preprocess=lambda data: process_data(data, ftype='pattern', dim=info['dim'], grid=info['grid']))
        field[varname] = field[varname] + data[varname].mean(dim='time')        


    # First EOF
    elif mode == 'first':
        
        filename = os.path.join(dirs['tmp'], f"{varname}_series_00000.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        timeseries = xr.open_mfdataset(filename, decode_times=time_coder, preprocess=lambda data: process_data(data, ftype='series', dim=info['dim'], grid=info['grid']))
        coeffs = timeseries.polyfit(dim='time', deg=1, skipna = True)
        theta = xr.polyval(xf, coeffs[f"{varname}_polyfit_coefficients"])        
        basis = pattern.isel(time=0)
        field += theta * basis
        field['time'] = xf

        # add mean (TO BE CHANGED!!!!!!!!!) USE THE READER
        logging.info(f"Adding mean trend of {varname}.")
        filename = os.path.join(dirs['tmp'], f"{varname}.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        data = xr.open_dataset(filename, decode_times=time_coder)
        #data = xr.open_mfdataset(filename) # use_cftime=True, preprocess=lambda data: process_data(data, ftype='pattern', dim=info['dim'], grid=info['grid']))
        field[varname] = field[varname] + data[varname].mean(dim='time')


    # Reconstruct the last frame (for dry-runs)
    elif mode == 'reco':

        for i in range(neofs):            
            filename = os.path.join(dirs['tmp'], f"{varname}_series_{str(i).zfill(5)}.nc")
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            timeseries = xr.open_mfdataset(filename, decode_times=time_coder, preprocess=lambda data: process_data(data, ftype='series', dim=info['dim'], grid=info['grid']))
            theta = timeseries[varname].isel(time=-1)
            basis = pattern.isel(time=i)
            field += theta * basis
        
        field = field.drop_vars({'time'})

        # add mean (TO BE CHANGED!!!!!!!!!)
        logging.info(f"Adding mean trend of {varname}.nc")
        filename = os.path.join(dirs['tmp'], f"{varname}.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        data = xr.open_dataset(filename, decode_times=time_coder)
        field[varname] = field[varname] + data[varname].mean(dim='time')  


    # Reconstruct the last frame (for dry-runs)
    elif mode == 'reco-first':
                    
        filename = os.path.join(dirs['tmp'], f"{varname}_series_00000.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        timeseries = xr.open_mfdataset(filename, decode_times=time_coder, preprocess=lambda data: process_data(data, ftype='series', dim=info['dim'], grid=info['grid']))
        theta = timeseries[varname].isel(time=-1)
        basis = pattern.isel(time=0)
        field += theta * basis

        field = field.drop_vars({'time'})

        # add mean (TO BE CHANGED!!!!!!!!!)
        logging.info(f"Adding mean trend of {varname}.nc")
        filename = os.path.join(dirs['tmp'], f"{varname}.nc")
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        data = xr.open_dataset(filename, decode_times=time_coder)
        field[varname] = field[varname] + data[varname].mean(dim='time')  


    # EOFs up to a percentage of the full
    elif mode == 'frac':
        
        threshold_percentage = 0.9

        weights = []
        cdo_command = f"cdo info -div {varname}_variance.nc -timsum {varname}_variance.nc"
        output = subprocess.check_output(cdo_command, shell=True, text=True)
        for line in output.splitlines():
            if "Mean" in line:
                mean_value = float(line.split(":")[3].strip())
                weights.append(mean_value)
        weights = np.array(weights)
        cumulative_weights = np.cumsum(weights) / np.sum(weights)

        # Find the index where cumulative sum exceeds the threshold percentage
        feofs = np.searchsorted(cumulative_weights, threshold_percentage) + 1

        for i in range(feofs):
            filename = os.path.join(dirs['tmp'], f"{varname}_series_{str(i).zfill(5)}.nc")
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            timeseries = xr.open_mfdataset(filename, decode_times=time_coder, preprocess=lambda data: process_data(data, ftype='series', dim=info['dim'], grid=info['grid']))
            p = timeseries.polyfit(dim='time', deg=1, skipna=True)
            theta = xr.polyval(xf, p[f"{varname}_polyfit_coefficients"])
            basis = pattern.isel(time=i)
            field += theta * basis

    # add linear regression fit modes: point-to-poit, global- & basin- based
    elif mode == 'fit':

        ds = reader_nemo_field(expname=expname, startyear=startyear, endyear=endyear, varname=varname)
        coeffs = ds.polyfit(dim='time', deg=1, skipna=True)
        field = xr.polyval(xf, coeffs[f"{varname}_polyfit_coefficients"])
        field['time'] = xf

        # add mean
        logging.info(f"Adding mean trend of {varname}.")
        field = field + ds.mean(dim='time')

    return field

##########################################################################################
# [8] forecast routines

def create_forecast_field(expname, varname, endyear, window, yearleap, 
                          mode='full', format='winter', smoothing=False, cleanup=False):
    """ 
    Function to forecast a single field using EOF
    
    Args:
    expname: experiment name
    varname: variable name
    endyear: final leg of the simulation
    window: years backward from endleg used by EOFs
    yearleap: years forward from endleg to forecast
    mode: EOF regression mode
    format: time format [plain, moving, winter, etc ...]
    smoothing: if needed, smooth out the forecasted fields
    
    """

    endleg = endyear - 1990 + 2
    startleg = endleg - window + 1
    startyear = 1990 + startleg - 2

    logging.info(f"Start/end year: {startyear}-{endyear}")
    logging.info(f"Time window: {window}")

    dirs = folders(expname)
    dirs['tmp'] = os.path.join(dirs['tmp'], f'Y{endyear}-W{window}-L{yearleap}')
    
    # prepare field and EOFs
    merge(expname, varname, startyear, endyear, window, yearleap)
    neofs=window-1    
    get_eofs(expname, varname, endyear, window, neofs, yearleap)

    # field projection in the future
    data = project_eofs(expname=expname, varname=varname, endyear=endyear, 
                        window=window, neofs=neofs, yearleap=yearleap, mode=mode)
    
    # apply constraints
    data = constraints_for_fields(data=data)

    # clean-up
    if cleanup:
        cleanup_eofs(expname=expname, varname=varname, endyear=endyear, window=window, yearleap=yearleap)

    return data


def forecaster(expname, varnames, endleg, window, yearleap, mode='full', format='winter', 
               smoothing=False, cleanup=False):
    """ 
    Function to assembly the forecast of multiple fields using EOF
    
    Args:
    expname: experiment name
    varname: variable name
    endleg: final leg of the simulation
    window: years backward from endleg used by EOFs
    yearleap: years forward from endleg to forecast
    mode: EOF regression mode
    smoothing: if needed, smooth out the forecasted fields
    
    """

    endyear = 1990 + endleg - 2

    # read forecast and change restart
    rdata = reader_rebuilt(expname, endleg, endleg)

    # create EOF
    for varname in varnames:
        
        field = create_forecast_field(expname, varname, endyear, window, yearleap, mode=mode, 
                                      format=format, smoothing=smoothing, cleanup=cleanup)

        field = field.rename({'time': 'time_counter', 'z': 'nav_lev'})
        field['time_counter'] = rdata['time_counter']

        # loop on the corresponding varlist    
        varlist = varlists.get(varname, []) # Get the corresponding varlist, default to an empty list if not found
        for vars in varlist:
            rdata[vars] = xr.where(rdata[vars] != 0.0, field[varname], 0.0)

    return rdata


def writer_restart(expname, rdata, leg):
    """ Write NEMO restart file in a temporary folder """

    dirs = folders(expname)
    flist = glob.glob(os.path.join(dirs['restart'], str(leg).zfill(3), expname + '*_' + 'restart' + '_????.nc'))
    timestep = get_nemo_timestep(flist[0])

    # ocean restart creation
    oceout = os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart.nc')
    rdata.to_netcdf(oceout, mode='w', unlimited_dims={'time_counter': True})

    # delete attributes
    delete_attrs(oceout)

    # copy ice restart
    orig = os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + timestep + '_restart_ice.nc')
    dest = os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart_ice.nc')
    shutil.copy(orig, dest)

    return None


def replacer(expname, leg):
    """ Replace modified restart files in the run folder """

    dirs = folders(expname)

    # cleaning
    browser = ['restart*.nc']
    for basefile in browser:
        filelist = sorted(glob.glob(os.path.join(dirs['exp'], basefile)))
        for file in filelist:
            if os.path.isfile(file):
                print('Removing' + file)
                os.remove(file)

    # create new links
    browser = ['restart.nc', 'restart_ice.nc']
    for file in browser:
        rebfile = os.path.join(dirs['tmp'], str(leg).zfill(3), file)
        resfile = os.path.join(dirs['restart'], str(leg).zfill(3), file)
        shutil.copy(rebfile, resfile)
        newfile = os.path.join(dirs['exp'], file)
        print("Linking rebuilt NEMO restart", file)            
        os.symlink(resfile, newfile)
    
    return None


def restorer(expname, leg):

    dirs = folders(expname)

    # copying from the restart folder required for the leg you asked
    browser = ['*restart*']
    for file in browser:
        filelist = sorted(glob.glob(os.path.join(dirs['restart'], str(leg).zfill(3), file)))
        for file in filelist:
            basefile = os.path.basename(file)
            targetfile = os.path.join(dirs['exp'], basefile)
            if not os.path.isfile(targetfile):
                if 'restart' in basefile:
                    newfile = os.path.join(dirs['exp'], '_'.join(basefile.split('_')[2:]))
                    print("Linking NEMO restart", file)
                    os.symlink(file, newfile)

    return None 

##########################################################################################
# [9] parser

def parse_args():
    """ Command line parser for nemo-restart """

    parser = argparse.ArgumentParser(description="Command Line Parser for NEMO forecast")

    # add positional argument (mandatory)
    parser.add_argument("expname", metavar="EXPNAME", help="Experiment name")
    parser.add_argument("leg", metavar="LEG", help="Leg of rebuilding", type=int)
    parser.add_argument("window", metavar="WINDOW", help="EOF window", type=int)
    parser.add_argument("yearleap", metavar="YEARLEAP", help="Year leap in the future", type=int)
    parser.add_argument("mode", metavar="MODE", help="EOF regression mode", type=str)

    # optional to activate nemo rebuild
    parser.add_argument("--rebuild", action="store_true", help="Enable NEMO rebuild")
    parser.add_argument("--forecast", action="store_true", help="Create NEMO forecast")
    parser.add_argument("--replace", action="store_true", help="Replace NEMO restarts")
    parser.add_argument("--restore", action="store_true", help="Restore nemo restart files")
    
    parsed = parser.parse_args()

    return parsed


if __name__ == "__main__":
    """
    Main entry point for the forecast_nemo script.

    This script processes the restart files for a specified experiment and time leg,
    creates necessary symbolic links, and runs the rebuild_nemo executable to rebuild
    the NEMO restart files.

    Command line arguments:
        expname (str): The name of the experiment.
        leg (str): The time leg of the experiment.
    """

    # parser
    args = parse_args()
    expname = args.expname
    leg = args.leg
    window = args.window
    yearleap = args.yearleap
    mode = args.mode

    # define folders
    dirs = folders(expname)

    # rebuild nemo restart files
    if args.rebuild:
        rebuilder(expname, leg)
    
    varnames = ['thetao', 'so']

    # forecast
    if args.forecast:
        rdata = forecaster(expname=expname, varnames=varnames, endleg=leg, window=window, yearleap=yearleap, mode=mode, smoothing=False)

        # write new restart file
        writer_restart(expname, rdata, leg)       

    # replace nemo restart files
    if args.replace:
        replacer(expname, leg)

    # restore
    if args.restore:
        restorer(expname, leg)

    # shutil.copy('logfile.log', os.path.join(dirs['tmp'], str(leg).zfill(3), 'logfile.log'))
