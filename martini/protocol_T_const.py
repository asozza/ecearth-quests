#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a command line tool to modfy the NEMO restart files from a specific EC-Eart4
experiment, given a specific experiment and leg. 

Needed modules:
# module load intel/2021.4.0 intel-mkl/19.0.5 prgenv/intel hdf5 netcdf4 
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/apps/netcdf4/4.9.1/INTEL/2021.4/lib:/usr/local/apps/hdf5/1.12.2/INTEL/2021.4/lib

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
import osprey_io as osi
import osprey_actions as osa
import osprey_means as osm


def parse_args():
    """Command line parser for nemo-restart"""

    parser = argparse.ArgumentParser(description="Command Line Parser for nemo-restart")

    # add positional argument (mandatory)
    parser.add_argument("expname", metavar="EXPNAME", help="Experiment name")
    parser.add_argument("leg", metavar="LEG", help="The leg you want to process for rebuilding", type=str)
    parser.add_argument("const", metavar="CONST", help="Constant delta temperature to be added (with sign)", type=float)

    # optional to activate nemo rebuild
    parser.add_argument("--rebuild", action="store_true", help="Enable nemo-rebuild")
    parser.add_argument("--replace", action="store_true", help="Replace nemo restart files")

    parsed = parser.parse_args()

    return parsed

if __name__ == "__main__":
    
    # parser
    args = parse_args()
    expname = args.expname
    leg = args.leg
    const = args.const

    # define folders
    dirs = osi.folders(expname)

    # rebuild nemo restart files
    if args.rebuild:
        osa.rebuild_nemo(expname, leg)

    # manipulate based on global temperature
    rdata = osa.manipulate_add_const(expname, leg, const)
    osi.write_nemo_restart(expname, rdata, leg)


