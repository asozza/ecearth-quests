#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a command line tool to modfy the NEMO restart files from a specific EC-Eart4
experiment, given a specific experiment and leg. 

Authors
Alessandro Sozza and Paolo Davini (CNR-ISAC, Nov 2023)
"""

import subprocess
import os
import glob
import shutil
import yaml
import argparse
import xarray as xr
from functions import preproc_nemo_T
from functions import moving_average
from functions import dateDecimal
from sklearn.linear_model import LinearRegression

# on atos, you will need to have
# module load intel/2021.4.0 intel-mkl/19.0.5 prgenv/intel hdf5 netcdf4 
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/apps/netcdf4/4.9.1/INTEL/2021.4/lib:/usr/local/apps/hdf5/1.12.2/INTEL/2021.4/lib

def parse_args():
    """Command line parser for nemo-restart"""

    parser = argparse.ArgumentParser(description="Command Line Parser for nemo-restart")

    # add positional argument (mandatory)
    parser.add_argument("expname", metavar="EXPNAME", help="Experiment name")
    parser.add_argument("leg", metavar="LEG", help="The leg you want to process for rebuilding", type=str)
    parser.add_argument("yearspan", metavar="YEARSPAN", help="Year span for fitting global temperature", type=int)
    parser.add_argument("yearleap", metavar="YEARLEAP", help="Year leap for projecting global temperature", type=int)

    # optional to activate nemo rebuild
    parser.add_argument("--rebuild", action="store_true", help="Enable nemo-rebuild")
    parser.add_argument("--replace", action="store_true", help="Replace nemo restart files")


    parsed = parser.parse_args()

    return parsed

def get_nemo_timestep(filename):
    """Minimal function to get the timestep from a nemo restart file"""

    return os.path.basename(filename).split('_')[1]

def rebuild_nemo(expname, leg, dirs):
    """Minimal nemo rebuilder in a temporary path"""

    rebuilder = os.path.join(dirs['rebuild'], "rebuild_nemo")
  
    for kind in ['restart', 'restart_ice']:
        print('Processing' + kind)
        flist = glob.glob(os.path.join(dirs['exp'], 'restart', leg.zfill(3), expname + '*_' + kind + '_????.nc'))
        tstep = get_nemo_timestep(flist[0])

        for filename in flist:
            destination_path = os.path.join(dirs['tmp'], os.path.basename(filename))
            try:
                os.symlink(filename, destination_path)
            except FileExistsError:
                pass

        rebuild_command = [rebuilder, os.path.join(dirs['tmp'],  expname + "_" + tstep + "_" + kind ), str(len(flist))]
        try:
            subprocess.run(rebuild_command, stderr=subprocess.PIPE, text=True, check=True)
            for file in glob.glob('nam_rebuld_*') : 
                os.remove(file)
        except subprocess.CalledProcessError as e:
            error_message = e.stderr
            print(error_message) 

        for filename in flist:
            destination_path = os.path.join(dirs['tmp'], os.path.basename(filename))
            os.remove(destination_path)


if __name__ == "__main__":
    
    # parser
    args = parse_args()
    expname = args.expname
    leg = args.leg
    yearspan = args.yearspan
    yearleap = args.yearleap

    # define directories
    dirs = {
        'exp': os.path.join("/ec/res4/scratch/itas/ece4", expname),
        'nemo': os.path.join("/ec/res4/scratch/itas/ece4/", expname, "output", "nemo"),
        'tmp':  os.path.join("/ec/res4/scratch/itas/martini", expname, leg.zfill(3)),
        'rebuild': "/ec/res4/hpcperm/itas/src/rebuild_nemo"
    }

    os.makedirs(dirs['tmp'], exist_ok=True)

    # rebuild nemo restart files
    if args.rebuild:
        rebuild_nemo(expname=expname, leg=leg, dirs=dirs)

    # extrapolate future global temperature
    legfile = os.path.join(dirs['exp'], 'leginfo.yml')
    with open(legfile, 'r', encoding='utf-8') as file:
        leginfo = yaml.load(file, Loader=yaml.FullLoader)
    info = leginfo['base.context']['experiment']['schedule']['leg']
    endyear = info['start'].year - 1
    startyear = endyear - yearspan
    print('Fitting in the range: ',startyear,endyear)

    # load domain
    domain = xr.open_dataset(os.path.join(dirs['exp'], 'domain_cfg.nc'))
    vol = domain['e1t']*domain['e2t']*domain['e3t_0']
    area = domain['e1t']*domain['e2t']

    # open database of the last chunk alone
    filelist = []
    for year in range(startyear, endyear):
        pattern = os.path.join(dirs['nemo'], f"{expname}_oce_1m_T_{year}-{year}.nc")
        matching_files = glob.glob(pattern)
        filelist.extend(matching_files)
    data = xr.open_mfdataset(filelist, preprocess=preproc_nemo_T)
    
    # extract averaged T filtered
    tt = dateDecimal(data['time'].values)
    toa = data['thetao'].weighted(vol).mean(dim=['z', 'y', 'x']).values.flatten()
    tom = moving_average(toa, 12)

    # fit global temperature
    Yg = [[tom[i]] for i in range(len(tom))]
    Xg = [[tt[i]] for i in range(len(tt))]
    model=LinearRegression()
    model.fit(Xg, Yg)
    mp = model.coef_[0][0]
    qp = model.intercept_[0]
    teq = mp*(endyear+yearleap) + qp
    print(' Fit coefficients: ',mp,qp)
    print(' Projected Temperature: ',teq)

    # modify restart files
    domain = domain.rename({'z': 'nav_lev'})
    vol = domain['e1t']*domain['e2t']*domain['e3t_0']
    filelist = glob.glob(os.path.join(dirs['tmp'],  expname + '*_restart.nc'))
    timestep = get_nemo_timestep(filelist[0])
    oce = os.path.join(dirs['tmp'], expname + '_' + timestep + '_restart.nc')
    xfield = xr.open_dataset(oce)
    varlist = ['tn', 'tb']
    for var in varlist:
        tef = xfield[var].where(xfield[var]!=0.0).isel(time_counter=0).weighted(vol).mean(dim=['nav_lev', 'y', 'x']).values
        print('Last value of Temperature: ',var,tef[0])
        xfield[var] = xr.where(xfield[var]!=0, xfield[var] - tef[0] + teq, 0.)
    
    # ocean restart creation
    oceout = os.path.join(dirs['tmp'], 'restart.nc')
    xfield.to_netcdf(oceout)

    # ice restart copy
    shutil.copy(os.path.join(dirs['tmp'], expname + '_' + timestep + '_restart_ice.nc'), os.path.join(dirs['tmp'], 'restart_ice.nc'))

    # replace nemo restart files
    if args.replace:
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
            rebfile = os.path.join(dirs['tmp'], file)
            resfile = os.path.join(dirs['exp'], 'restart', leg.zfill(3), file)
            shutil.copy(rebfile, resfile)
            newfile = os.path.join(dirs['exp'], file)
            print("Linking rebuilt NEMO restart", file)            
            os.symlink(resfile, newfile)

