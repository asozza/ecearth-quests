#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reader module

Author: Alessandro Sozza
Date: March 2024
"""

import os
import glob
import shutil
import logging
import numpy as np
import xarray as xr

from osprey.means.means import timemean, spacemean
from osprey.utils.utils import get_nemo_timestep
from osprey.utils.folders import folders
from osprey.utils.time import get_leg, dateDecimal
from osprey.actions.rebuilder import rebuilder


##########################################################################################
# Readers of NEMO output

def _nemodict(grid, freq):
    """Dictionary of NEMO output fields"""

    gridlist = ["T", "U", "V", "W"]    
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
    """General preprocessing routine for NEMO data based on grid type"""
    
    grid_mappings = _nemodict(grid, None)[grid]  # None for freq as it is not used here

    data = data.rename_dims({grid_mappings["x_grid"]: 'x', grid_mappings["y_grid"]: 'y'})
    data = data.rename({
        grid_mappings["nav_lat"]: 'lat', 
        grid_mappings["nav_lon"]: 'lon', 
        grid_mappings["depth"]: 'z', 
        'time_counter': 'time'
    })
    data = data.swap_dims({grid_mappings["x_grid_inner"]: 'x', grid_mappings["y_grid_inner"]: 'y'})
    data = data.drop_vars(['time_centered'], errors='ignore')
    data = data.drop_dims(['axis_nbounds'], errors='ignore')

    return data


def preproc_nemo_ice(data):
    """Preprocessing routine for NEMO for ice"""

    data = data.rename({'time_counter': 'time'})
    
    return data


def read_nemo(expname, startyear, endyear, grid="T", freq="1m"):
    """Main function to read nemo data"""

    dirs = folders(expname)
    dict = _nemodict(grid, freq)

    filelist = []
    for year in range(startyear, endyear+1):
        pattern = os.path.join(dirs['nemo'], f"{expname}_{dict[grid]['format']}_{year}-{year}.nc")
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    logging.info('Files to be loaded %s', filelist)
    data = xr.open_mfdataset(filelist, preprocess=lambda d: dict[grid]["preproc"](d, grid), use_cftime=True)

    return data


##########################################################################################
# Reader of NEMO domain

def preproc_nemo_domain(data):
    """ preprocessing routine for nemo domain """

    data = data.rename({'time_counter': 'time'})

    return data

def read_domain(expname):
    """ read NEMO domain configuration file """

    dirs = folders(expname)
    filename = os.path.join(dirs['exp'], 'domain_cfg.nc')
    domain = xr.open_mfdataset(filename, preprocess=preproc_nemo_domain)
    domain = domain.isel(time=0)

    return domain

def elements(expname):
    """ define differential forms for integrals """

    df = {}
    domain = read_domain(expname=expname)
    df['vol'] = domain['e1t']*domain['e2t']*domain['e3t_0']
    df['area'] = domain['e1t']*domain['e2t']
    df['dx'] = domain['e1t']
    df['dy'] = domain['e2t']
    df['dz'] = domain['e3t_0']

    return df


##########################################################################################
##########################################################################################
# Readers of averaged data: timeseries, profiles, maps, hovmoller, pdfs etc ...

def read_timeseries_T(expname, startyear, endyear, var):
    """ read averaged timeseries (for generic grid) """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"timeseries_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_profile_T(expname, startyear, endyear, var):
    """ read averaged profiles """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"profiles_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_hovmoller_T(expname, startyear, endyear, var):
    """ read averaged hovmoller diagram """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"hovmoller_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_map_T(expname, startyear, endyear, var):
    """read averaged map """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"map_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_field_T(expname, startyear, endyear, var):
    """read averaged meanfield in 2D or 3D """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"field_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

#### for anomaly

def read_timeseries_local_anomaly_T(expname, startyear, endyear, var):
    """ read averaged rms local anomaly timeseries """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"timeseries_local_anomaly_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_profile_local_anomaly_T(expname, startyear, endyear, var):
    """ read averaged profiles of rms local anomaly """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"profiles_local_anomaly_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_hovmoller_local_anomaly_T(expname, startyear, endyear, var):
    """ read averaged rms local anomaly hovmoller diagram """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"hovmoller_local_anomaly_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data

def read_map_local_anomaly_T(expname, startyear, endyear, var):
    """read averaged map of rms local anomaly """

    dirs = folders(expname)
    filename = os.path.join(dirs['perm'], f"map_local_anomaly_{var}_{startyear}-{endyear}.nc")
    data = xr.open_dataset(filename, use_cftime=True)

    return data


##########################################################################################
# Containers for reader/creator of averaged data

# for averaged timeseries
def read_averaged_timeseries_T(expname, startyear, endyear, inivar, ndim, isub, iload):
    """ reader/creator of averaged timeseries for T grid """

    # check if var is a new variable defined on a subregion
    if '-' in inivar:
        var = inivar.split('-')[0]
    else:
        var = inivar

    dirs = folders(expname)
    df = elements(expname)

    # try to read averaged data
    try:
        data = read_timeseries_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    #if iload == 'orig':
    data = read_T(expname, startyear, endyear)
    print(" Averaging ... ")
    tvec = dateDecimal(data['time'].values)
    vec = spacemean(expname, data[var], ndim)
    #elif iload == 'post':
    #    data = read_from_cdo_T(expname, startyear, endyear, var)
    #    tvec = ost.dateDecimal(data['time'].values)
    #    vec = osm.spacemean(expname, data[var], ndim)

    # ONLY for 3D variables! 
    # compute var in subregions: mixed layer (mix), pycnocline (pyc), abyss (aby)
    if isub == True:
        subvec = spacemean3d_suball(expname, data[var])
        ds = xr.Dataset({
            'time': xr.DataArray(data = tvec, dims = ['time'], coords = {'time': tvec}, 
                            attrs = {'units' : 'years', 'long_name' : 'years'}),
            var : xr.DataArray(data = vec, dims = ['time'], coords = {'time': tvec}, 
                            attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name}),                            
            var+'-mix' : xr.DataArray(data = subvec[0], dims = ['time'], coords = {'time': tvec}, 
                            attrs  = {'units' : data[var].units, 'long_name' : 'mixed layer ' + data[var].long_name}), 
            var+'-pyc' : xr.DataArray(data = subvec[1], dims = ['time'], coords = {'time': tvec}, 
                            attrs  = {'units' : data[var].units, 'long_name' : 'pycnocline ' + data[var].long_name}), 
            var+'-aby' : xr.DataArray(data = subvec[2], dims = ['time'], coords = {'time': tvec}, 
                            attrs  = {'units' : data[var].units, 'long_name' : 'abyss ' + data[var].long_name})},
            attrs = {'description': 'ECE4/NEMO 1D timeseries averaged from T_grid variables'})
        
        # write the averaged data and read it again
        print(" Saving averaged data ... ")
        filename = os.path.join(dirs['perm'], f"timeseries_{var}_{startyear}-{endyear}.nc")
        ds.to_netcdf(filename)
        data = read_timeseries_T(expname, startyear, endyear, var)

        return data

    # create xarray dataset
    ds = xr.Dataset({
        'time': xr.DataArray(data = tvec, dims = ['time'], coords = {'time': tvec}, 
                        attrs = {'units' : 'years', 'long_name' : 'years'}), 
        var : xr.DataArray(data = vec, dims = ['time'], coords = {'time': tvec}, 
                        attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})},
        attrs = {'description': 'ECE4/NEMO 1D timeseries averaged from T_grid variables'})

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"timeseries_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_timeseries_T(expname, startyear, endyear, var)

    return data

# for averaged profiles
def read_averaged_profile_T(expname, startyear, endyear, var, iload):
    """ reader/creator of averaged profiles for T grid """

    dirs = folders(expname)
    df = elements(expname)

    # try to read averaged data
    try:
        data = read_profile_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    #if iload == 'orig':
    data = read_T(expname, startyear, endyear)
    print(" Averaging ... ")
    # and spatial averaging of the desidered variable
    zvec = data['z'].values.flatten()
    vec = data[var].weighted(df['area']).mean(dim=['time', 'y', 'x']).values.flatten()
    #elif iload == 'post':
    #    data = read_from_cdo_T(expname, startyear, endyear, var)
    #    zvec = data['z'].values.flatten()
    #    vec = data[var].weighted(df['area']).mean(dim=['time', 'y', 'x']).values.flatten()

    # create xarray dataset
    ds = xr.Dataset({
        'z': xr.DataArray(data = zvec, dims = ['z'], coords = {'z': zvec}, 
                             attrs = {'units' : 'm', 'long_name' : 'depth'}), 
        var : xr.DataArray(data = vec, dims = ['z'], coords = {'z': zvec}, 
                               attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})}, 
        attrs = {'description': 'ECE4/NEMO 1D averaged profiles from T_grid variables'})

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"profiles_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)    
    data = read_profile_T(expname, startyear, endyear, var)

    return data

# for averaged hovmöller diagram
def read_averaged_hovmoller_T(expname, startyear, endyear, var):
    """ reader/creator of averaged hovmöller diagram for T grid """

    dirs = folders(expname)
    df = elements(expname)

    # try to read averaged data
    try:
        data = read_hovmoller_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    data = read_T(expname, startyear, endyear)
    print(" Averaging ... ")
    # and spatial averaging of the desidered variable
    tvec = dateDecimal(data['time'].values)
    vec = spacemean(expname, data[var], '2D')

    # create xarray dataset
    print(" Allocating new xarray dataset ... ")
    ds = xr.Dataset({
        'time': xr.DataArray(data = tvec, dims = ['time'], coords = {'time': tvec}, 
                    attrs = {'units' : 'years', 'long_name' : 'years'}), 
        'z': xr.DataArray(data = data['z'], dims = ['z'], coords = {'z': data['z']}, 
                    attrs = {'units' : 'm', 'long_name' : 'depth'}), 
        var : xr.DataArray(data = vec, dims = ['time', 'z'], coords = {'time': tvec, 'z': data['z']},
                        attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})}, 
            attrs = {'description': 'ECE4/NEMO averaged T_grid 2D hovmoller diagram'})

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"hovmoller_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_field_T(expname, startyear, endyear, var)

    return data

# for averaged field
def read_averaged_field_T(expname, startyear, endyear, var, ndim):
    """ reader/creator of averaged field for T grid """

    dirs = folders(expname)
    df = elements(expname)

    # try to read averaged data
    try:
        data = read_field_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    data = read_T(expname, startyear, endyear)
    print(" Averaging ... ")
    # and spatial averaging of the desidered variable
    vec = timemean(data[var])

    # create xarray dataset
    print(" Allocating new xarray dataset ... ")
    if ndim == '3D':
        ds = xr.Dataset({
            'lat': xr.DataArray(data = data['lat'], dims = ['y', 'x'], coords = {'y': data['y'], 'x': data['x']}, 
                        attrs = {'units' : 'deg', 'long_name' : 'latitude'}),
            'lon': xr.DataArray(data = data['lon'], dims = ['y', 'x'], coords = {'y': data['y'], 'x': data['x']}, 
                        attrs = {'units' : 'deg', 'long_name' : 'longitude'}),                   
            'z': xr.DataArray(data = data['z'], dims = ['z'], coords = {'z': data['z']},
                        attrs = {'units' : 'm', 'long_name' : 'depth'}),                             
            var : xr.DataArray(data = vec, dims = ['z', 'y', 'x'], coords = {'z': data['z'], 'y': data['y'], 'x': data['x']},
                        attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})}, 
            attrs = {'description': 'ECE4/NEMO averaged T_grid_3D field'})

    if ndim == '2D':
        ds = xr.Dataset({
            'lat': xr.DataArray(data = data['lat'], dims = ['y', 'x'], coords = {'y': data['y'], 'x': data['x']}, 
                        attrs = {'units' : 'deg', 'long_name' : 'latitude'}),
            'lon': xr.DataArray(data = data['lon'], dims = ['y', 'x'], coords = {'y': data['y'], 'x': data['x']}, 
                        attrs = {'units' : 'deg', 'long_name' : 'longitude'}),                   
            var : xr.DataArray(data = vec, dims = ['y', 'x'], coords = {'y': data['y'], 'x': data['x']},
                        attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})}, 
            attrs = {'description': 'ECE4/NEMO averaged T_grid_2D field'})

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"field_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_field_T(expname, startyear, endyear, var)

    return data


##################################################################################################################
# Reader/Creator of averaged local anomaly

# timeseries of averaged rms local anomaly
def read_averaged_timeseries_local_anomaly_T(expname, startyear, endyear, refname, startref, endref, var, ndim):
    """ reader/creator of averaged timeseries of rms local anomaly for T grid """

    dirs = folders(expname)
    df = elements(expname)

    # try to read averaged data
    try:
        data = read_timeseries_local_anomaly_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    # read field and meanfield
    mdata = read_averaged_field_T(refname, startref, endref, var, ndim)
    data = read_T(expname, startyear, endyear)
    delta = np.square(data[var]-mdata[var]) # choose cost function: square
    
    print(" Averaging ... ")
    # spatial averaging of the desidered variable
    tvec = dateDecimal(data['time'].values)
    vec = spacemean(expname, delta, ndim)
    vec = np.power(vec,0.5)

    # create xarray dataset
    ds = xr.Dataset({
        'time': xr.DataArray(data = tvec, dims = ['time'], coords = {'time': tvec}, 
                        attrs = {'units' : 'years', 'long_name' : 'years'}), 
        var : xr.DataArray(data = vec, dims = ['time'], coords = {'time': tvec}, 
                        attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})},
        attrs = {'description': 'ECE4/NEMO 1D timeseries averaged rms local anomaly from T_grid variables'})

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"timeseries_local_anomaly_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_timeseries_local_anomaly_T(expname, startyear, endyear, var)

    return data

# profile of the anomaly 
def read_averaged_profile_local_anomaly_T(expname, startyear, endyear, refname, startref, endref, var):
    """ reader/creator of averaged profile of rms local anomaly for T grid """

    dirs = folders(expname)
    df = elements(expname)    

    # try to read averaged data
    try:
        data = read_profile_local_anomaly_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    # read field and meanfield
    mdata = read_averaged_field_T(refname, startref, endref, var, '3D')
    data = read_T(expname, startyear, endyear)
    delta = np.power(data[var]-mdata[var],2) # choose cost function: square

    print(" Averaging ... ")
    # and spatial averaging of the desidered variable
    zvec = data['z'].values.flatten()
    vec = delta.weighted(df['area']).mean(dim=['time', 'y', 'x']).values.flatten()
    vec = np.power(vec,0.5)

    # create xarray dataset
    print(" Allocating new xarray dataset ... ")
    ds = xr.Dataset({
        'z': xr.DataArray(data = data['z'], dims = ['z'], coords = {'z': data['z']}, 
                    attrs = {'units' : 'm', 'long_name' : 'depth'}), 
        var : xr.DataArray(data = vec, dims = ['z'], coords = {'z': data['z']},
                    attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})
                    }, 
        attrs = {'description': 'ECE4/NEMO averaged T_grid 2D profile of rms local anomaly'}
        )

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"hovmoller_local_anomaly_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_hovmoller_local_anomaly_T(expname, startyear, endyear, var)

    return data

# howmoller diagram of averaged rms local anomaly
def read_averaged_hovmoller_local_anomaly_T(expname, startyear, endyear, refname, startref, endref, var):
    """ reader/creator of averaged hovmöller diagram of rms local anomaly for T grid """

    dirs = folders(expname)
    
    # try to read averaged data
    try:
        data = read_hovmoller_local_anomaly_T(expname, startyear, endyear, var)
        print(" Averaged data found ")
        return data
    except FileNotFoundError:
        print(" Averaged data not found. Creating new file ... ")

    # If averaged data not existing, read original data
    print(" Loading data ... ")
    # read field and meanfield
    mdata = read_averaged_field_T(refname, startref, endref, var, '3D')
    data = read_T(expname, startyear, endyear)
    delta = np.power(data[var]-mdata[var],2) # choose cost function: square

    print(" Averaging ... ")
    # and spatial averaging of the desidered variable
    tvec = dateDecimal(data['time'].values)
    vec = spacemean(expname, delta, '2D')
    vec = np.power(vec,0.5)

    # create xarray dataset
    print(" Allocating new xarray dataset ... ")
    ds = xr.Dataset({
        'time': xr.DataArray(data = tvec, dims = ['time'], coords = {'time': tvec}, 
                    attrs = {'units' : 'years', 'long_name' : 'years'}), 
        'z': xr.DataArray(data = data['z'], dims = ['z'], coords = {'z': data['z']}, 
                    attrs = {'units' : 'm', 'long_name' : 'depth'}), 
        var : xr.DataArray(data = vec, dims = ['time', 'z'], coords = {'time': tvec, 'z': data['z']},
                    attrs  = {'units' : data[var].units, 'long_name' : data[var].long_name})
                    }, 
        attrs = {'description': 'ECE4/NEMO averaged T_grid 2D hovmoller diagram of rms local anomaly'}
        )

    # write the averaged data and read it again
    print(" Saving averaged data ... ")
    filename = os.path.join(dirs['perm'], f"hovmoller_local_anomaly_{var}_{startyear}-{endyear}.nc")
    ds.to_netcdf(filename)
    data = read_hovmoller_local_anomaly_T(expname, startyear, endyear, var)

    return data


##########################################################################################
# Reader/Writer of NEMO restart files

def read_restart(expname, startyear, endyear):
    """ Reader of NEMO restart files in a range of legs """

    startleg = get_leg(startyear)
    endleg = get_leg(endyear)
    dirs = folders(expname)

    try:
        data = read_rebuilt(expname, startleg, endleg)
        return data
    except FileNotFoundError:
        print(" Restart file not found. Rebuilding ... ")

    # rebuild files
    for leg in range(startleg,endleg+1):
        rebuilder(expname, leg)

    data = read_rebuilt(expname, startleg, endleg)

    return data

def read_rebuilt(expname, startleg, endleg):
    """ Read rebuilt NEMO restart files """

    dirs = folders(expname)
    
    filelist = []
    for leg in range(startleg,endleg+1):
        pattern = os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '*_restart.nc')
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    data = xr.open_mfdataset(filelist, use_cftime=True)

    return data

def write_restart(expname, rdata, leg):
    """ Write NEMO restart file """

    dirs = folders(expname)
    flist = glob.glob(os.path.join(dirs['restart'], str(leg).zfill(3), expname + '*_' + 'restart' + '_????.nc'))
    timestep = get_nemo_timestep(flist[0])

    # ocean restart creation
    oceout = os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart.nc')
    rdata.to_netcdf(oceout, mode='w', unlimited_dims={'time_counter': True})

    # copy ice restart
    orig = os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + timestep + '_restart_ice.nc')
    dest = os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart_ice.nc')
    shutil.copy(orig, dest)

    return None

##########################################################################################

