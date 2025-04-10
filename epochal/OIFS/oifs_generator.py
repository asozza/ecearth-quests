#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a command line tool to OIFS ICs and BCs from default available ones.
It can produce data from using CDO and GRIB_API. 
It uses cdo bindings for python in a rough way to allow for exploration of temporary files.


Authors
Paolo Davini (CNR-ISAC, Apr 2024)
"""

import subprocess
import os
import shutil
import cdo
from utils import extract_grid_info, ecmwf_grid
cdo = cdo.Cdo()
cdo.debug = True

# configurable
target_grid = 'TL63L31'
BASE_TGT = '/lus/h2resw01/scratch/ccpd/TL63-cy48'

do_clean = False
do_spectral = True
do_boundary = True
do_gridpoint = True

#-----------------------#
# configurable with caution
startdate = '19900101'
source_grid = 'TL255L91'
# Higher resolution is better, in principle.
# However, given the issue that we have with remapcon, we will use coarser resolution for now.

# where original OIFS data is found
OIFS_BASE = '/home/ccpd/hpcperm/ECE4-DATA/oifs/'

# these are only available on ATOS
OIFS_BC='/perm/smw/oifs48-icmcl'

# climate files version (climate.v020 for cy48, v015 for cy43)
climate_version = "climate.v020"

# these are available on atos and can be produced by oifs_create_corners.py
GRIDS = '/home/ccpd/perm/ecearth4/oifs-grids'

# temporart directory
TMPDIR = '/ec/res4/scratch/ccpd/tmpic_cy48_check'


#---------------------------#

def generate_spectral_conditions(OIFS_IC, spectral, TMPDIR):
    """
    Generate spectral initial conditions.
    """
    # This is done with a clean spectral truncation with cdo.
    # Orography is therefore realiable.
    # The file has to be split in two since orography is GRIB1 and the rest is GRIB2
    print("Truncating spectral file ICMSHECE4INIT to", spectral, "harmonics")
    grib2file = cdo.sp2sp(spectral, input=f"-selname,lnsp,vo,t,d {OIFS_IC}/ICMSHECE4INIT")
    gribtemp = cdo.selname("z", input=f"{OIFS_IC}/ICMSHECE4INIT", options="--eccodes")
    grib1file = cdo.sp2sp(spectral, input=gribtemp, options="--eccodes")
    subprocess.call(f"cat {grib2file} {grib1file} > {TMPDIR}/ICMSHECE4INIT", shell=True)

def generate_gridpoint_conditions(OIFS_IC, GRIDS, target_spectral, TMPDIR):
    """
    Generate gridpoint initial conditions.
    """
    print("Remapping ICMGG gaussian ICs to target grid")
    for file in ["ICMGGECE4INIT", "ICMGGECE4INIUA"]:
        # This is the alternative to be done with remapcon using the grid fils computed with oifs_create_corner.py
        #icmtmp = cdo.remapcon(f"{GRIDS}/{target_spectral}_grid.nc",
        #             input=f"-setgrid,{GRIDS}/{source_spectral}_grid.nc {OIFS_IC}/{file}")
        #cdo.setgrid(f"grids/{target_spectral}.txt", input=icmtmp, output=f"{TMPDIR}/{file}"
        # this is the old version with remapnn
        cdo.remapnn(f"{GRIDS}/{target_spectral}_grid.nc", input=f"{OIFS_IC}/{file}", output=f"{TMPDIR}/{file}")

def generate_boundary_conditions(OIFS_BC, target_grid, climate_version, TMPDIR, BC_TGT, ecmwf_name):
    """
    Generate boundary conditions.
    """
    climfile = f"{OIFS_BC}/{target_grid}/{climate_version}/ICMCLECE4"
    # case when BCs are already available, produced by Klaus Wyser
    if os.path.exists(climfile):
        print("BCs already exist")
        shutil.copy(climfile, f"{BC_TGT}/ICMCLECE4")
    # This is done with a mergetime of the 7 variables in the ECMWF directory based on a magic command by Klaus Wyser
    else:
        print("Building BCs from", ecmwf_name, "data")
        variables = ["alb", "aluvp", "aluvd", "alnip", "alnid", "lail", "laih"]
        paths = [f"{OIFS_BC}/{ecmwf_name}/month_{var}" for var in variables]
        cdo.mergetime(options="-L", input=paths, output=f"{TMPDIR}/temp.grb")
        cdo.settaxis("2021-01-15,00:00:00,1month",
                     input=f"{TMPDIR}/temp.grb",
                     output=f"{BC_TGT}/ICMCLECE4-1990")
        os.remove(f"{TMPDIR}/temp.grb")

def vertical_interpolation(TMPDIR, GRIDS, target_spectral, vertical, IC_TGT, do_clean):
    """
    Perform vertical interpolation for spectral and gridpoint data.
    """

    # Procedure for vertical interpolation requires all the data to be in grid point space.
    # This is done by converting the spectral fields to gaussian grids and then moving back them to the spectral space 
    # It has been decided to interpolate spectral data (T, D, V) and keep gaussian data (Q, etc.) on the gaussian reduced grid
    # Orography and surface pressure are not touched and attached to the files at the end of the operations
    # A-B coefficients for remapeta are downloaded from ECMWF website and then converted to txt file 
    # in CDO-compliant style with convert_aka_bika.py script. These are stored in the grids folder. 
    # To set gaussian reduced grids the grid files are produced with descriptor_generator.py and
    # also stored in txt file in the grids folder
    print("Vertical interpolation is necessary")
    print("Select z and lnsp from ICMSHECE4INIT...")
    VERTVALUES = (int(vertical) + 1) * 2
    cdo.selname("z", input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/orog.grb")
    cdo.selname("lnsp", input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/lnsp.grb")
    
    # Modify lnsp to avoid CDO issues
    subprocess.call(f"grib_set -s numberOfVerticalCoordinateValues={VERTVALUES} {TMPDIR}/lnsp.grb {TMPDIR}/lnsp2.grb", shell=True)

    # Convert spectral fields to Gaussian grid
    print("Converting ICMSHECE4INIT to Gaussian grid")
    cdo.sp2gpl(input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/sp2gauss.grb")

    # Remap spectral fields to Gaussian reduced grid
    print("Remapping spectral fields from Gaussian to Gaussian reduced")
    remapped = cdo.remapcon(f"{GRIDS}/{target_spectral}_grid.nc", input=f"{TMPDIR}/sp2gauss.grb")
    cdo.setgrid(f"grids/{target_spectral}.txt", input=remapped, output=f"{TMPDIR}/sp2gauss_reduced.grb")

    # Merge files for interpolation
    print("Merging files")
    subprocess.call(f"cat {TMPDIR}/ICMGGECE4INIUA {TMPDIR}/sp2gauss_reduced.grb > {TMPDIR}/single.grb", shell=True)

    # Perform hybrid levels interpolation
    print("Remapping vertical on hybrid levels")
    gridfile = f"grids/L{vertical}.txt"
    cdo.remapeta(gridfile, input=f"{TMPDIR}/single.grb", output=f"{TMPDIR}/remapped.grb")

    # Create INITUA file
    print("Selecting fields to create ICMSHECE4INIUA and setting Gaussian reduced grid")
    cdo.setgrid(f"grids/{target_spectral}.txt", input=f"-selname,q,o3,crwc,cswc,clwc,ciwc,cc {TMPDIR}/remapped.grb", output=f"{IC_TGT}/ICMGGECE4INIUA")

    # Convert back to spectral space
    print("Converting back to spectral (through Gaussian regular) the spectral fields")
    cdo.gp2spl(input=f"-setgridtype,regular -selname,t,vo,d {TMPDIR}/remapped.grb", output=f"{TMPDIR}/spback.grb")

    # Merge with orography and lnsp to create the final SH file
    print("Merging files and creating the final ICMSHECE4INIT")
    subprocess.call(f"cat {TMPDIR}/spback.grb {TMPDIR}/lnsp2.grb {TMPDIR}/orog.grb > {IC_TGT}/ICMSHECE4INIT", shell=True)

    # Move the ICMGGECE4INIT file
    shutil.move(f"{TMPDIR}/ICMGGECE4INIT", f"{IC_TGT}/ICMGGECE4INIT")

    # Cleanup temporary files if required
    if do_clean:
        print("Cleaning up")
        for file in ["gp2gauss.grb", "sp2gauss.grb", "sp2gauss_reduced.grb", "single.grb", "remapped.grb",
                     "ICMSHECE4INIT", "ICMGGECE4INIUA", "spback.grb", "lnsp.grb", "lnsp2.grb", "orog.grb"]:
            os.remove(f"{TMPDIR}/{file}")

if __name__ == "__main__":

    IC_TGT = os.path.join(BASE_TGT, target_grid, startdate)
    BC_TGT = os.path.join(BASE_TGT, target_grid, 'climate')

    for d in [IC_TGT, BC_TGT, TMPDIR]:
        os.makedirs(d, exist_ok=True)

    # target grid info
    grid_type, spectral, vertical = extract_grid_info(target_grid)
    ecmwf_name = str(spectral) + ecmwf_grid(grid_type)
    target_spectral = 'T' + grid_type + str(spectral)

    # source grid info
    ic_grid_type, ic_spectral, ic_vertical = extract_grid_info(source_grid)
    source_spectral = 'T' + ic_grid_type + str(ic_spectral)

    OIFS_IC = os.path.join(OIFS_BASE, source_grid, startdate)

    if ic_vertical != vertical:
        print("Vertical interpolation is necessary")
        DO_VERTICAL = True
    else:
        DO_VERTICAL = False

    if do_spectral:
        generate_spectral_conditions(OIFS_IC, spectral, TMPDIR)

    # Gridpoint generation
    if do_gridpoint:
        OIFS_IC = os.path.join(OIFS_BASE, source_grid, startdate)
        generate_gridpoint_conditions(OIFS_IC, GRIDS, target_spectral, TMPDIR)

    # Boundary condition generation
    if do_boundary:
        generate_boundary_conditions(OIFS_BC, target_grid, climate_version, TMPDIR, BC_TGT, ecmwf_name)


    # move the files to the target directory
    if do_spectral and do_boundary:
        if not DO_VERTICAL:
            print("Copying files to the target directory")
            for file in ["ICMSHECE4INIT", "ICMGGECE4INIT", "ICMGGECE4INIUA"]:
                shutil.move(f"{TMPDIR}/{file}", f"{IC_TGT}/{file}")


        else:
            vertical_interpolation(TMPDIR, GRIDS, target_spectral, vertical, IC_TGT, do_clean)

    if do_clean:
        os.rmdir(TMPDIR)

    print("Done")
