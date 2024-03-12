#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 _____  __________  ____  ______   __ 
/  __ \/   __   _ \|  _ \|  ___ \ / / 
| |  |    |_ | |_|   |_|   |__ \ V /  
| |  | |\_  \|  __/|    /|  __| | |   
| |__| |__|    |   | |\ \| |____| |   
\___________/|_|   |_| \__________|   

OSPREY: Ocean Spin-uP acceleratoR for Earth climatologY
--------------------------------------------------------
Osprey library for i/o operations

Authors
Alessandro Sozza (CNR-ISAC, Mar 2024)
"""

import subprocess
import os
import glob
import yaml
import numpy as np
import xarray as xr
import cftime

def folders(expname, leg):

    dirs = {
        'exp': os.path.join("/ec/res4/scratch/itas/ece4", expname),
        'nemo': os.path.join("/ec/res4/scratch/itas/ece4/", expname, "output", "nemo"),
        'tmp':  os.path.join("/ec/res4/scratch/itas/martini", expname, leg.zfill(3)),
        'rebuild': "/ec/res4/hpcperm/itas/src/rebuild_nemo",
        'backup': os.path.join("/ec/res4/scratch/itas/ece4", expname + "-backup")
    }
    os.makedirs(dirs['tmp'], exist_ok=True)

    return dirs

def preproc_nemo_T(data):
    """preprocessing routine for nemo for T grid"""

    data = data.rename_dims({'x_grid_T': 'x', 'y_grid_T': 'y'})
    data = data.rename({'deptht': 'z', 'time_counter': 'time', 'thetao': 'to'})
    data = data.swap_dims({'x_grid_T_inner': 'x', 'y_grid_T_inner': 'y'})
    data.coords['z'] = -data['z']
    
    return data

def preproc_nemo_U(data):
    """preprocessing routine for nemo for U grid"""

    data = data.rename_dims({'x_grid_U': 'x', 'y_grid_U': 'y'})
    data = data.rename({'depthu': 'z', 'time_counter': 'time'})
    data = data.swap_dims({'x_grid_U_inner': 'x', 'y_grid_U_inner': 'y'})
    data.coords['z'] = -data['z']
    
    return data

def preproc_nemo_V(data):
    """preprocessing routine for nemo for V grid"""

    data = data.rename_dims({'x_grid_V': 'x', 'y_grid_V': 'y'})
    data = data.rename({'depthv': 'z', 'time_counter': 'time'})
    data = data.swap_dims({'x_grid_V_inner': 'x', 'y_grid_V_inner': 'y'})
    data.coords['z'] = -data['z']
    
    return data

def preproc_nemo_W(data):
    """preprocessing routine for nemo for W grid"""

    data = data.swap_dims({'x_grid_W_inner': 'x', 'y_grid_W_inner': 'y'})
    data = data.rename_dims({'x_grid_W': 'x', 'y_grid_W': 'y'})
    data = data.rename({'depthw': 'z', 'time_counter': 'time'})
    data.coords['z'] = -data['z']

    return data

def preproc_nemo_ice(data):
    """preprocessing routine for nemo for ice"""

    data = data.rename({'time_counter': 'time'})
    
    return data

def readmf_T(expname, startyear, endyear):

    dirs = folders(expname)
    filelist = []
    for year in range(startyear, endyear):
        pattern = os.path.join(dirs['nemo'], f"{expname}_oce_*_T_{year}-{year}.nc")
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    data = xr.open_mfdataset(filelist, preprocess=preproc_nemo_T, use_cftime=True)

    return data

def readmf_ice(expname, startyear, endyear):

    dirs = folders(expname)
    filelist = []
    for year in range(startyear, endyear):
        pattern = os.path.join(dirs['nemo'], f"{expname}_ice_*_{year}-{year}.nc")
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    data = xr.open_mfdataset(filelist, preprocess=preproc_nemo_ice, use_cftime=True)

    return data

def read_T(expname, year):

    dirs = folders(expname)
    filelist = os.path.join(dirs['nemo'], f"{expname}_oce_*_T_{year}-{year}.nc")
    data = xr.open_mfdataset(filelist, preprocess=preproc_nemo_T, use_cftime=True)

    return data

def read_ice(expname, year):

    dirs = folders(expname)
    filelist = os.path.join(dirs['nemo'], f"{expname}_ice_*_{year}-{year}.nc")
    data = xr.open_mfdataset(filelist, preprocess=preproc_nemo_ice, use_cftime=True)

    return data

def read_domain(expname):

    dirs = folders(expname)
    domain = xr.open_dataset(os.path.join(dirs['exp'], 'domain_cfg.nc'))

    return domain

def get_nemo_timestep(filename):
    """Minimal function to get the timestep from a nemo restart file"""

    return os.path.basename(filename).split('_')[1]

