#!/usr/bin/env python3
"""
Process NEMO PALEORCA2 bathymetry data by opening and closing specific regions.

This script takes an input directory containing bathymetry NetCDF files,
applies regional modifications, and saves the processed output.
"""

import argparse
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt


def close_region(xfield, region):
    """Close specific regions by setting bathymetry to 0."""
    if region == 'Caspian':
        xfield['bathy_metry'][:, 129:140,  154:159] = 0
    if region == 'Victoria':
        xfield['bathy_metry'][:,  95:102,  148] = 0
    if region == "GreatLakes":
        xfield['bathy_metry'][:,  132:141,  86:94] = 0
    if region == "Arctic":
        xfield['bathy_metry'][:,  173,  117] = 0
        xfield['bathy_metry'][:,  173,  63:65] = 0
        xfield['bathy_metry'][:,  172,  103:106] = 0
        #xfield['bathy_metry'][:,  173,  110:116] = 0
    if region == "Britain":
        xfield['bathy_metry'][:,  140:142,  130] = 0
    if region == "Panama":
        xfield['bathy_metry'][:,  119:120,  84] = 0
        xfield['bathy_metry'][:,  115,  89:93] = 0
        xfield['bathy_metry'][:,  114,  93] = 0
    if region == "Thailand":
        # xfield['bathy_metry'][:,  105:108, 1] = 0
        # xfield['bathy_metry'][:,  103, 2] = 0
        # xfield['bathy_metry'][:,  114:116, 1] = 0
        # xfield['bathy_metry'][:,  110:114, 2] = 0
        # xfield['bathy_metry'][:,  107:110, 3] = 0
        xfield['bathy_metry'][:,  112:113, 2] = 0
    if region == "Italy":
        xfield['bathy_metry'][:,  133:134,  138] = 0
        xfield['bathy_metry'][:,  132,  139] = 0
    if region == "Kamchatka":
        xfield['bathy_metry'][:,  146:149,  35] = 0
        xfield['bathy_metry'][:,  148,  36] = 0
    if region == "Barents":
        xfield['bathy_metry'][:,  153:154, 144:147] = 0
    if region == "RedSea":
        xfield['bathy_metry'][:,  123,  149] = 0.
    return xfield


def open_region(xfield, region):
    """Open specific regions by setting bathymetry to specific depths."""
    if region == "Gibraltair":
        xfield['bathy_metry'][:,  130,  131] = 284.
    if region == "RedSea":
        xfield['bathy_metry'][:,  116,  153:155] = 137.
    if region == "Hormuz":
        xfield['bathy_metry'][:,  122,  159] = 100.
    if region == "Adriatic":
        xfield['bathy_metry'][:,  132,  140] = 300.
    if region == "Baltic":
        xfield['bathy_metry'][:,  149,  140] = 50.
    if region == "Kara":
        xfield['bathy_metry'][:,  165,  149] = 25.
        xfield['bathy_metry'][:,  165,  152] = 25.
    if region == "Australia":
        xfield['bathy_metry'][:,  83:84,  23] = 100.
    if region == "Greenland":
        xfield['bathy_metry'][:,  168,  121] = 300.
    if region == "Black":
        xfield['bathy_metry'][:,  133,  145] = 200.
    if region == "Indonesia":
        xfield['bathy_metry'][:,  102,  11] = 200.
    
    return xfield


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process NEMO bathymetry data by opening and closing specific regions"
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="Input directory containing bathymetry NetCDF files"
    )
    parser.add_argument(
        "--infile",
        type=str,
        default="PALEORCA_bathy_metry_from_eORCA1_nn.nc",
        help="Input bathymetry file name (default: PALEORCA_bathy_metry_from_eORCA1_nn.nc)"
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default="PALEORCA_bathy_metry_from_eORCA1_nn_closed.nc",
        help="Output bathymetry file name (default: PALEORCA_bathy_metry_from_eORCA1_nn_closed.nc)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate a plot of the processed bathymetry"
    )
    
    return parser.parse_args()


def main():
    """Main processing function."""
    args = parse_arguments()
    
    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist.")
        sys.exit(1)
    
    # Construct file paths
    input_file = os.path.join(args.input_dir, args.infile)
    output_file = os.path.join(args.input_dir, args.outfile)
    
    # Validate input file
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        sys.exit(1)
    
    print(f"Processing bathymetry data from: {input_file}")
    print(f"Output will be saved to: {output_file}")
    
    try:
        # Load the dataset
        print("Loading bathymetry data...")
        xfield = xr.open_dataset(input_file)
        
        # Apply regional closures (from the notebook workflow)
        print("Applying regional closures...")
        xfield = close_region(xfield, 'Caspian')
        xfield = close_region(xfield, 'Victoria')
        xfield = close_region(xfield, 'GreatLakes')
        xfield = close_region(xfield, 'RedSea')
        xfield = close_region(xfield, 'Thailand')
        
        # Apply regional openings (from the notebook workflow)
        print("Applying regional openings...")
        xfield = open_region(xfield, 'Gibraltair')
        xfield = open_region(xfield, 'Hormuz')
        xfield = open_region(xfield, 'Baltic')
        xfield = open_region(xfield, 'Adriatic')
        xfield = open_region(xfield, 'Kara')
        xfield = open_region(xfield, 'Australia')
        xfield = open_region(xfield, 'Black')
        xfield = open_region(xfield, 'Greenland')
        xfield = open_region(xfield, 'Indonesia')
        
        # Re-open Red Sea as done in the notebook
        xfield = open_region(xfield, 'RedSea')
        
        # Generate plot if requested
        if args.plot:
            print("Generating plot...")
            plt.figure(figsize=(12, 8))
            xfield['bathy_metry'].plot(vmax=1)
            plt.title("Processed Bathymetry Data")
            plot_file = os.path.join(args.input_dir, "bathymetry_plot.png")
            plt.savefig(plot_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {plot_file}")
            plt.close()
        
        # Save the processed dataset
        print("Saving processed data...")
        xfield.to_netcdf(output_file)
        
        print("Processing complete!")
        print(f"Processed bathymetry saved to: {output_file}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()