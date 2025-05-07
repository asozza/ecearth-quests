#!/usr/bin/env python3

import os
import shutil
import argparse
from datetime import datetime
import subprocess
import yaml

BASEDIR='/home/ccpd/scratch/ece4'

def parse_arguments():
    parser = argparse.ArgumentParser(description="Rollback script for resetting experiment data.")
    parser.add_argument("expname", type=str, help="Name of the experiment (e.g., 'ABT4').")
    parser.add_argument("--years", type=int, default=5, help="Number of years to rollback.")
    parser.add_argument("--clean", action="store_true", help="Enable cleaning and restoring from backup.")
    parser.add_argument("--backup", action="store_true", help="Create a backup of the files before modifying them.")
    parser.add_argument("--dry", action="store_true", help="Print commands without executing them.")
    return parser.parse_args()

def read_current_date_from_yaml(yaml_file):
    """Reads the current date from the YAML file."""
    with open(yaml_file, "r") as file:
        data = yaml.safe_load(file)
    # Extract the date from the YAML structure
    current_date = data["base.context"]["experiment"]["schedule"]["leg"]["start"]
    return current_date
    #return datetime.strptime(current_date, "%Y-%m-%d %H:%M:%S").strftime("%Y%m%d")

def compute_adjusted_date(current_date, rollback_years):
    """
    Computes the adjusted date
    """
    # Convert the dates to datetime objects
    yyyy = int(datetime.strftime(current_date, "%Y"))
    mm = datetime.strftime(current_date, "%m")

    # Calculate the difference
    newyear = yyyy - rollback_years
    adjusted_date = datetime.strptime(f'{newyear}{mm}01', "%Y%m%d")
    #adjusted_date = current_date_dt - delta + timedelta(days=1)
    # Return the adjusted date in YYYYMMDD format
    return adjusted_date

def create_backup(file_path):
    """Creates a backup of the specified file."""
    backup_path = f"{file_path}.backup"
    if os.path.exists(file_path):
        print(f"Creating backup: {backup_path}")
        shutil.copy(file_path, backup_path)

def clean_and_restore(file_path):
    """Cleans and restores a file from its backup."""
    backup_path = f"{file_path}.backup"
    if os.path.exists(backup_path):
        print(f"Restoring {file_path} from backup...")
        shutil.copy(backup_path, file_path)
    else:
        print(f"No backup found for {file_path}. Skipping...")

def rollback(expname, clean, backup, dry, rollback_years):
    DIR = os.path.join(BASEDIR, expname)
    yaml_file = os.path.join(DIR, "leginfo.yml")
    current_date = read_current_date_from_yaml(yaml_file)
    print(f"Current date: {current_date}")

    #original_date = subprocess.run(['cdo', '-s', 'showdate', f"{DIR}/ICMGG{expname}INIUA"], capture_output=True, text=True)
    #original_date = original_date.stdout.strip()
    #original_date = datetime.strptime(original_date, "%Y-%m-%d")
    #print(f"Original date: {original_date}")

    adjusted_date = compute_adjusted_date(current_date, rollback_years=rollback_years)
    print(f"Adjusted date: {adjusted_date}")

    if adjusted_date > current_date:
        print(f"Rollback date {adjusted_date} is later than current date {current_date}.")
        return

    # Compute the difference in seconds and days
    diff_sec = int((current_date - adjusted_date).total_seconds())
    diff_days = (current_date - adjusted_date).days

    CSTEP = diff_sec // 3600
    CTIME = f"{diff_days:08d}0000"
    print(f"CSTEP: {CSTEP} ; CTIME={CTIME}")

    exit()

    grib_filenames = [
        f"{DIR}/ICMGG{expname}INIUA",
        f"{DIR}/ICMSH{expname}INIT",
        f"{DIR}/ICMGG{expname}INIT"
    ]

    if backup:
        print("Creating backups...")
        for grib_file in grib_filenames:
            create_backup(grib_file)
        create_backup(f"{DIR}/rcf")
        for file in os.listdir(DIR):
            if file.startswith("srf"):
                create_backup(file)
        
    if clean:
        print("Cleaning and restoring from backup...")
        for grib_file in grib_filenames:
            clean_and_restore(grib_file)
        clean_and_restore(f"{DIR}/rcf")
        for file in os.listdir(DIR):
            if file.startswith("srf"):
                clean_and_restore(os.path.join(DIR, file))

   
    print("Changing the date of the initial conditions...")

    for filename in grib_filenames:
        old_file = f"{filename}.old"
        print(f"grib_set -s dataDate={adjusted_date} {old_file} {filename}")
        if os.path.exists(filename):
            if not dry:
                os.rename(filename, old_file)
                subprocess.run(["grib_set", "-s", f"dataDate={adjusted_date}", old_file, filename])

    # Modify OIFS restart control file
    print("Modifying the old rcf file...")

    rcf_file = f"{DIR}/rcf"
    if os.path.exists(rcf_file):
        with open(rcf_file, "r") as file:
            rcf_content = file.read()

        rcf_content = rcf_content.replace(
            "CSTEP[^,]*,", f'CSTEP   = "    {CSTEP}",'
        )
        rcf_content = rcf_content.replace(
            "CTIME[^,]*,", f'CTIME   = "{CTIME:<9}  ",'
        )

        if not dry:
            with open(rcf_file, "w") as file:
                file.write(rcf_content)

    print(f"CSTEP: {CSTEP} ; CTIME={CTIME}")

    print("Renaming the restart files...")
    for file in os.listdir(DIR):
        if file.startswith("srf"):
            ext = file.split(".")[-1]
            new_file = f"srf{CTIME}.{ext}"
            print(f"cp {file} {new_file}")
            if not dry:
                shutil.copy(file, new_file)
                os.remove(file)

if __name__ == "__main__":
    args = parse_arguments()
    rollback(args.expname, args.clean, args.backup, args.dry, args.years)