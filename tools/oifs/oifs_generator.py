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
import argparse
import shutil
import cdo
from utils import extract_grid_info, ecmwf_grid
from utils import unpack_grib_file, repack_grib_file
cdo = cdo.Cdo()
cdo.debug = True

# configurable
target_grid = 'TL63L31'
BASE_TGT = '/home/ccpd/hpcperm/ECE4-DATA/oifs/'

do_clean = False
do_spectral = True
do_boundary = False
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
# climate files version (climate.v020 for cy48, v015 for cy43)
climate_version = "climate.v020"
OIFS_BC = {
    "ready": "/perm/smw/oifs48-icmcl",
    "compute": os.path.join('/home/rdx/data/climate/', climate_version)
}


# these are available on atos and can be produced by oifs_create_corners.py
GRIDS = '/home/ccpd/perm/ecearth4/oifs-grids'

# temporart directory
TMPDIR = '/ec/res4/scratch/ccpd/tmpic_cy48'

# CDO grib2 options
GRIB2="-f grb2 --eccodes"
GRIB1="-f grb1 --eccodes"


#---------------------------#

def parse_arguments():
    """
    Parse command-line arguments for the script.
    """
    parser = argparse.ArgumentParser(description="Generate OIFS ICs and BCs with configurable target resolution.")
    parser.add_argument(
        "-t", "--target", type=str, required=True,
        help="Target grid resolution (e.g., TL63L31, TL255L91)."
    )
    parser.add_argument(
        "--clean", action="store_true",
        help="Clean up temporary files after processing."
    )

    return parser.parse_args()

def generate_spectral_conditions(OIFS_IC, spectral, TMPDIR):
    """
    Generate spectral initial conditions.
    This is done with a clean spectral truncation with cdo.
    Orography is reliable and should not produce Gibbs oscillations.
    """
    # NEW PROCEDURE FOR CY48
    # in cycle all is grb2 so no problem with this
    print("Truncating spectral file ICMSHECE4INIT to", spectral, "harmonics")
    cdo.sp2sp(spectral, input=f"{OIFS_IC}/ICMSHECE4INIT", output=f"{TMPDIR}/ICMSHECE4INIT", options=GRIB2)
    
    # OLD PROCEDURE FOR CY43
    # The file has to be split in two since orography is GRIB1 and the rest is GRIB2
    #grib2file = cdo.sp2sp(spectral, input=f"-selname,lnsp,vo,t,d {OIFS_IC}/ICMSHECE4INIT")
    #gribtemp = cdo.selname("z", input=f"{OIFS_IC}/ICMSHECE4INIT", options="--eccodes")
    #grib1file = cdo.sp2sp(spectral, input=gribtemp, options="--eccodes")
    #subprocess.call(f"cat {grib2file} {grib1file} > {TMPDIR}/ICMSHECE4INIT", shell=True)
    #if clean:
    #    os.remove(grib2file)
    #    os.remove(grib1file)


def generate_gridpoint_conditions(OIFS_IC, target_spectral, TMPDIR, clean=True):
    """
    Generate gridpoint initial conditions.
    Nearest neighbour interpolation is used to remap the data.
    """
    print("Remapping ICMGG gaussian ICs to target grid")
    for file in ["ICMGGECE4INIT", "ICMGGECE4INIUA"]:
        # This is the alternative to be done with remapcon using the grid fils computed with oifs_create_corner.py
        #icmtmp = cdo.remapcon(f"{GRIDS}/{target_spectral}_grid.nc",
        #             input=f"-setgrid,{GRIDS}/{source_spectral}_grid.nc {OIFS_IC}/{file}")
        #cdo.setgrid(f"grids/{target_spectral}.txt", input=icmtmp, output=f"{TMPDIR}/{file}")
        # this is the old version with remapnn
        #cdo.remapnn(f"{GRIDS}/{target_spectral}_grid.nc", input=f"{OIFS_IC}/{file}", output=f"{TMPDIR}/{file}", options="--eccodes")
        gridfile = f"grids/{target_spectral}.txt"
        inputfile = f"{OIFS_IC}/{file}"
        # Split GRIB messages in GRIB1 and GRIB2 to avoid CDO problems in cy48
        grib1, grib2 = unpack_grib_file(inputfile, f'{TMPDIR}/{file}')

        # interpolate
        for grbfile in [grib1, grib2]:
            if os.path.exists(grbfile):
                cdo.remapnn(gridfile, 
                            input=grbfile,
                            output=f"{grbfile}.remap", options="--eccodes")   
        # pack them together
        repack_grib_file(f"{grib1}.remap", f"{grib2}.remap", f"{TMPDIR}/{file}")
        if clean:
            for file in [grib1, grib2]:
                if os.path.exists(file):
                    os.remove(file)
                    os.remove(f"{file}.remap")

def generate_boundary_conditions(OIFS_BC, target_grid, climate_version, TMPDIR, BC_TGT, ecmwf_name):
    """
    Generate boundary conditions.
    """
    climfile = f"{OIFS_BC['ready']}/{target_grid}/{climate_version}/ICMCLECE4"
    # case when BCs are already available, produced by Klaus Wyser
    if os.path.exists(climfile):
        print("BCs already exist")
        shutil.copy(climfile, f"{BC_TGT}/ICMCLECE4")
    # This is done with a mergetime of the 7 variables in the ECMWF directory based on a magic command by Klaus Wyser
    else:
        print("Building BCs from", ecmwf_name, "data")
        variables = ["alb", "aluvp", "aluvd", "alnip", "alnid", "lail", "laih"]
        paths = [f"{OIFS_BC['compute']}/{ecmwf_name}/month_{var}" for var in variables]
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
    cdo.selname("z", input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/orog.grb", options=GRIB2)
    cdo.selname("lnsp", input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/lnsp.grb", options=GRIB2)
    
    # Modify lnsp to avoid CDO issues
    subprocess.call(f"grib_set -s numberOfVerticalCoordinateValues={VERTVALUES} {TMPDIR}/lnsp.grb {TMPDIR}/lnsp2.grb", shell=True)

    # Convert spectral fields to Gaussian grid
    print("Converting ICMSHECE4INIT to Gaussian grid")
    cdo.sp2gpl(input=f"{TMPDIR}/ICMSHECE4INIT", output=f"{TMPDIR}/sp2gauss.grb", options=GRIB2)

    # Remap spectral fields to Gaussian reduced grid
    print("Remapping spectral fields from Gaussian to Gaussian reduced")
    remapped = cdo.remapcon(f"{GRIDS}/{target_spectral}_grid.nc", input=f"{TMPDIR}/sp2gauss.grb", options=GRIB2)
    cdo.setgrid(f"grids/{target_spectral}.txt", input=remapped, output=f"{TMPDIR}/sp2gauss_reduced.grb")

    # Merge files for interpolation
    print("Merging files")
    subprocess.call(f"grib_copy {TMPDIR}/ICMGGECE4INIUA {TMPDIR}/sp2gauss_reduced.grb {TMPDIR}/single.grb", shell=True)
    #subprocess.call(f"cat {TMPDIR}/ICMGGECE4INIUA {TMPDIR}/sp2gauss_reduced.grb > {TMPDIR}/single.grb", shell=True)

    # Perform hybrid levels interpolation
    print("Remapping vertical on hybrid levels")
    gridfile = f"grids/L{vertical}.txt"
    cdo.remapeta(gridfile, input=f"{TMPDIR}/single.grb", output=f"{TMPDIR}/remapped.grb", options=GRIB2)

    # Create INITUA file
    print("Selecting fields to create ICMSHECE4INIUA and setting Gaussian reduced grid")
    cdo.setgrid(f"grids/{target_spectral}.txt", input=f"-selname,q,o3,crwc,cswc,clwc,ciwc,cc {TMPDIR}/remapped.grb", output=f"{IC_TGT}/ICMGGECE4INIUA", options=GRIB2)

    # Convert back to spectral space
    print("Converting back to spectral (through Gaussian regular) the spectral fields")
    cdo.gp2spl(input=f"-setgridtype,regular -selname,t,vo,d {TMPDIR}/remapped.grb", output=f"{TMPDIR}/spback.grb", options="--eccodes")

    # Merge with orography and lnsp to create the final SH file
    print("Merging files and creating the final ICMSHECE4INIT")
    #subprocess.call(f"cat {TMPDIR}/spback.grb {TMPDIR}/lnsp2.grb {TMPDIR}/orog.grb > {IC_TGT}/ICMSHECE4INIT", shell=True)
    subprocess.call(f"grib_copy {TMPDIR}/spback.grb {TMPDIR}/lnsp2.grb {TMPDIR}/orog.grb {IC_TGT}/ICMSHECE4INIT", shell=True)

    # Move the ICMGGECE4INIT file
    shutil.move(f"{TMPDIR}/ICMGGECE4INIT", f"{IC_TGT}/ICMGGECE4INIT")

    # Cleanup temporary files if required
    if do_clean:
        print("Cleaning up")
        for file in ["sp2gauss.grb", "sp2gauss_reduced.grb", "single.grb", "remapped.grb",
                     "ICMSHECE4INIT", "ICMGGECE4INIUA", "spback.grb", "lnsp.grb", "lnsp2.grb", "orog.grb"]:
            os.remove(f"{TMPDIR}/{file}")

if __name__ == "__main__":

    # Parse command-line arguments
    args = parse_arguments()
    target_grid = args.target
    do_clean = args.clean

    IC_TGT = os.path.join(BASE_TGT, target_grid, startdate)
    BC_TGT = os.path.join(BASE_TGT, target_grid, climate_version)

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
        generate_gridpoint_conditions(OIFS_IC, target_spectral, TMPDIR, clean=do_clean)

    # Boundary condition generation
    if do_boundary:
        generate_boundary_conditions(OIFS_BC, target_grid, climate_version, TMPDIR, BC_TGT, ecmwf_name)


    # move the files to the target directory
    if do_spectral or do_boundary:
        if not DO_VERTICAL:
            print("Copying files to the target directory")
            for file in ["ICMSHECE4INIT", "ICMGGECE4INIT", "ICMGGECE4INIUA"]:
                shutil.move(f"{TMPDIR}/{file}", f"{IC_TGT}/{file}")


        else:
            vertical_interpolation(TMPDIR, GRIDS, target_spectral, vertical, IC_TGT, do_clean)

    if do_clean:
        os.rmdir(TMPDIR)

    print("Done")
