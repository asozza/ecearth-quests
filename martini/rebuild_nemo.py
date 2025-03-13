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
import yaml
import shutil
import platform
import argparse
import subprocess


def load_config():
    """ Load configuration file """

    config_path = "../../config.yaml"
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    return {}
    
def folders(expname):
    """ List of global paths dependent on expname """
    
    system = platform.system().lower()
    config = load_config()
    config = config.get(system)    
    base_path = config.get("base_path")
    src_path = config.get("src_path")
    data_path = config.get("data_path")

    dirs = {
        'exp': os.path.join(base_path, expname),
        'nemo': os.path.join(base_path, expname, "output", "nemo"),
        'oifs': os.path.join(base_path, expname, "output", "oifs"),
        'restart': os.path.join(base_path, expname, "restart"),
        'log': os.path.join(base_path, expname, "log"),
        'tmp': os.path.join(base_path, expname, "tmp"),
        'post': os.path.join(base_path, expname, "post"),
        'rebuild': os.path.join(src_path, "rebuild_nemo"),
        'domain': os.path.join(data_path, "nemo", "domain")
    }

    # Create 'post' & 'tmp' folder if it doesn't exist
    if not os.path.exists(dirs['post']):
        os.makedirs(dirs['post'])

    if not os.path.exists(dirs['tmp']):
        os.makedirs(dirs['tmp'])

    return dirs

def get_nemo_timestep(filename):
    """ Get timestep from a NEMO restart file """

    return os.path.basename(filename).split('_')[1]


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

    rebuild_exe = os.path.join(dirs['rebuild'], "rebuild_nemo.sh")
  
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
            print(rebuild_command)
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
    #tstep = get_nemo_timestep(glob.glob(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '*_restart.nc'))[0])
    shutil.copy(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + tstep + '_restart.nc'), os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart.nc'))
    shutil.copy(os.path.join(dirs['tmp'], str(leg).zfill(3), expname + '_' + tstep + '_restart_ice.nc'), os.path.join(dirs['tmp'], str(leg).zfill(3), 'restart_ice.nc'))

    # delete temporary files
    flist = glob.glob('nam_rebuild*')
    for file in flist:
        os.remove(file)

    return None


def parse_args():
    """Command line parser for rebuild_nemo

    This function parses command line arguments for the rebuild_nemo script.
    It expects two positional arguments: the experiment name and the time leg.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """

    parser = argparse.ArgumentParser(description="Command Line Parser for rebuild_nemo")

    # add positional argument (mandatory)
    parser.add_argument("expname", metavar="EXPNAME", help="Experiment name")
    parser.add_argument("leg", metavar="LEG", help="Time leg", type=str)
    parsed = parser.parse_args()

    return parsed

if __name__ == "__main__":
    """
    Main entry point for the rebuild_nemo script.

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

    rebuilder(expname, leg)
