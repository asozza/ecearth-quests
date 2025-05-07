#!/usr/bin/env python3

import os
import shutil
import argparse
from datetime import datetime
import subprocess
import yaml

BASEDIR='/home/ccpd/scratch/ece4'
ROLLBACK_DATE = '20001001'  # YYYYMMDD format

def parse_arguments():
    parser = argparse.ArgumentParser(description="Rollback script for resetting experiment data.")
    parser.add_argument("expname", type=str, help="Name of the experiment (e.g., 'ABT4').")
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
    print(f"Current date from YAML: {current_date}")
    return current_date.strftime("%Y%m%d")
    #return datetime.strptime(current_date, "%Y-%m-%d %H:%M:%S").strftime("%Y%m%d")

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

def rollback(expname, clean, backup, dry):
    DIR = os.path.join(BASEDIR, expname)
    yaml_file = os.path.join(DIR, "leginfo.yml")
    current_date = read_current_date_from_yaml(yaml_file)
    

    # Convert dates to datetime objects
    date1 = datetime.strptime(ROLLBACK_DATE, "%Y%m%d")
    date2 = datetime.strptime(current_date, "%Y%m%d")

    if date1 > date2:
        print(f"Rollback date {ROLLBACK_DATE} is later than current date {current_date}.")
        return

    # Compute the difference in seconds and days
    diff_sec = int((date2 - date1).total_seconds())
    diff_days = (date2 - date1).days

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
        print(f"grib_set -s dataDate={ROLLBACK_DATE} {old_file} {filename}")
        if os.path.exists(filename):
            if not dry:
                os.rename(filename, old_file)
                subprocess.run(["grib_set", "-s", f"dataDate={ROLLBACK_DATE}", old_file, filename])

    # Modify OIFS restart control file
    print("Modifying the old rcf file...")
    CSTEP = diff_sec // 3600
    CTIME = f"{diff_days:08d}0000"

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
    rollback(args.expname, args.clean, args.backup, args.dry)