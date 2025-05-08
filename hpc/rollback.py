#!/usr/bin/env python3

import os
import re
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
    parser.add_argument("--restore", action="store_true", help="Full restore of the experiment.")
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
    """Creates a backup of the specified file or symlink."""
    backup_path = f"{file_path}.backup"
    if os.path.exists(file_path):
        if os.path.islink(file_path):
            # Se il file è un symlink, crea un symlink di backup
            target = os.readlink(file_path)  # Ottieni il percorso a cui punta il symlink
            print(f"Creating symlink backup: {backup_path} -> {target}")
            if not os.path.exists(backup_path):
                os.symlink(target, backup_path)  # Crea un nuovo symlink per il backup
        else:
            # Se il file non è un symlink, copia normalmente
            print(f"Creating backup: {backup_path}")
            shutil.copy(file_path, backup_path, follow_symlinks=True)

def clean_and_restore(file_path):
    """Cleans and restores a file or symlink from its backup."""
    backup_path = f"{file_path}.backup"
    if os.path.exists(backup_path):
        if os.path.islink(backup_path):
            # Se il backup è un symlink, ripristina il symlink
            target = os.readlink(backup_path)  # Ottieni il percorso a cui punta il symlink
            print(f"Restoring symlink: {file_path} -> {target}")
            if os.path.exists(file_path):
                os.remove(file_path)  # Rimuovi il file originale
            os.symlink(target, file_path)  # Ripristina il symlink
        else:
            # Se il backup non è un symlink, copia normalmente
            print(f"Restoring file: {file_path} from {backup_path}")
            shutil.copy(backup_path, file_path)
    else:
        print(f"No backup found for {file_path}. Skipping...")

def rollback(expname, restore, backup, dry, rollback_years):
    DIR = os.path.join(BASEDIR, expname)
    
    backup_files = [file for file in os.listdir(DIR) if file.endswith(".backup")]
    if not backup_files and not backup:
        print("No backup files found. Please create backups first.")
        return
    
    if restore and backup:
        print("Restore and backup options cannot be used together.")
        return
    
    if backup_files and backup:
        print("Backup files already exist. Please remove them before creating new backups.")
        return

    yaml_file = os.path.join(DIR, "leginfo.yml")
    current_date = read_current_date_from_yaml(yaml_file)
    print(f"Current date: {current_date}")

    original_date = subprocess.run(['cdo', '-s', 'showdate', f"{DIR}/ICMGG{expname}INIUA"], capture_output=True, text=True)
    original_date = original_date.stdout.strip()
    original_date = datetime.strptime(original_date, "%Y-%m-%d")
    print(f"Original date: {original_date}")

    adjusted_date = compute_adjusted_date(current_date, rollback_years=rollback_years)
    print(f"Adjusted date: {adjusted_date}")

    if adjusted_date > current_date:
        print(f"Rollback date {adjusted_date} is later than current date {current_date}.")
        return

    # Compute the difference in seconds and days
    #diff_sec = (current_date-original_date).total_seconds() - int((current_date - adjusted_date).total_seconds())
    #diff_days = (current_date-original_date).days - (current_date - adjusted_date).days
    diff_sec = int((current_date - adjusted_date).total_seconds())
    diff_days = (current_date - adjusted_date).days

    CSTEP = diff_sec // 3600
    CTIME = f"{diff_days:08d}0000"
    CSTOP = f"d{(diff_days+365):07d}"
    NRESTS = CSTEP + 8760
    print(f"CSTEP: {CSTEP} ; CTIME={CTIME}")
    print(f"CSTOP: {CSTOP} ; NRESTS={NRESTS}")
    grib_filenames = [
        f"{DIR}/ICMGG{expname}INIUA",
        f"{DIR}/ICMSH{expname}INIT",
        f"{DIR}/ICMGG{expname}INIT"
    ]
    rcf = f"{DIR}/rcf"
    namelist = os.path.join(DIR, "templates", "namelist.oifs.j2")


    if backup:
        print("Creating backups...")
        for grib_file in grib_filenames:
            create_backup(grib_file)
        create_backup(rcf)
        for file in os.listdir(DIR):
            if file.startswith("srf"):
                file_path = os.path.join(DIR, file) 
                create_backup(file_path)
        create_backup(namelist)
        print("Backups created successfully.")
        return

        
    print("Cleaning and restoring from backup...")
    for grib_file in grib_filenames:
        clean_and_restore(grib_file)
    clean_and_restore(rcf)
    clean_and_restore(namelist)
    for file in os.listdir(DIR):
        if file.startswith("srf") and file.endswith(".backup"):
            file_path = os.path.join(DIR, file)
            original_path = file_path.replace(".backup", "")
            print(f"Cleaning and restoring {file_path} to {original_path}")
            if os.path.exists(original_path):
                os.remove(original_path)
            os.symlink(os.readlink(file_path), original_path)
    
    if restore:
        for file in os.listdir(DIR):
            if file.endswith(".backup"):
                print(f"Removing backup file: {file}")
                os.remove(os.path.join(DIR, file))
            if CTIME in file:
                print(f"Removing file: {file}")
                os.remove(os.path.join(DIR, file))
        print("Restoration completed successfully.")
        return

    print("Changing the date of the initial conditions...")
    adjusted_date_str = adjusted_date.strftime("%Y%m%d")
    for filename in grib_filenames:
        old_file = f"{filename}.old"
        print(f"grib_set -s dataDate={adjusted_date_str} {old_file} {filename}")
        if os.path.exists(filename):
            if not dry:
                os.rename(filename, old_file)
                subprocess.run(["grib_set", "-s", f"dataDate={adjusted_date_str}", old_file, filename])
                os.remove(old_file)
    # Modify OIFS restart control file
    print("Modifying the old rcf file...")
    
    if os.path.exists(rcf):
        with open(rcf, "r") as file:
            rcf_content = file.read()

        # Sostituzione per CSTEP
        rcf_content = re.sub(
            r"CSTEP[^,]*,", 
            f'CSTEP   = "    {CSTEP}",', 
            rcf_content
        )

        # Sostituzione per CTIME
        rcf_content = re.sub(
            r"CTIME[^,]*,", 
            f'CTIME   = "{CTIME:<9}  ",', 
            rcf_content
        )
        print(rcf_content)
        if not dry:
            with open(rcf, "w") as file:
                file.write(rcf_content)

    print("Modyfing the namelist template...")
    # Pattern to find the jinja-style expression and replace it with 'd0000730'
    cstop_pattern = r"d\{\{\s*\"%07d\"\s*%\s*\(experiment\.schedule\.leg\.end\s*-\s*experiment\.schedule\.start\)\.days\s*\}\}"
    nrests_pattern = r"\{\{\s*\(24\*3600\*\(experiment\.schedule\.leg\.end\s*-\s*experiment\.schedule\.start\)\.days\s*/\s*model_config\.oifs\.grid\.dt\)\s*\|\s*int\s*\}\}"

    with open(namelist, 'r') as file:
        content = file.read()
    content = re.sub(cstop_pattern, CSTOP, content)
    content = re.sub(nrests_pattern, str(NRESTS), content)
    with open(namelist, 'w') as file:
        file.write(content)

    print("Renaming the restart files...")
    for file in os.listdir(DIR):
        if CTIME in file:
            print(f"Removing {DIR}/{file}...")
            os.remove(os.path.join(DIR, file))
        if file.startswith("srf") and not CTIME in file and not file.endswith(".backup"):
            ext = file.split(".")[-1]
            new_file = f"{DIR}/srf{CTIME}.{ext}"
            print(f"Linking restart {DIR}/{file} {new_file}...")
            if not dry:
                if os.path.exists(new_file):
                    print(f"Removing existing file: {new_file}")
                    os.remove(new_file)
                os.symlink(os.readlink(os.path.join(DIR,file)), new_file)
                os.remove(os.path.join(DIR,file))

    print("Rollback completed. Finger crossed for the next run.")

if __name__ == "__main__":
    args = parse_arguments()
    rollback(args.expname, args.restore, args.backup, args.dry, args.years)
