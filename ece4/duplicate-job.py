#!/usr/bin/env python3
"""
Simple tool to duplicate a EC-Earth job. 
This script duplicates an existing job by copying its directory and renaming the configuration files.

Usage:
    python duplicate-job.py -i/--expname1 <source_expname> -t/--expname2 <target_expname> -c/--config <config_file>

    -i/--expname1: Name of the source experiment (e.g., aa00).
    -o/--expname2: Name of the target experiment (e.g., bb00).
    -c/--config: Path to the configuration file (default: config.yml).
"""

import os
import shutil
import argparse

from yaml_util import load_yaml

def duplicate_job(expname1, expname2, config, clean=False):
    """
    Duplicate an existing job to create a new one with a different experiment ID.

    Args:
        expname1 (str): Name of the source experiment.
        expname2 (str): Name of the target experiment.
        config (str): Path to the configuration file.
    """
    conf = load_yaml(config, expand_env=True)
    job_dir = conf['job_dir']

    source_dir = os.path.join(job_dir, expname1)
    target_dir = os.path.join(job_dir, expname2)

    if not os.path.exists(source_dir):
        raise ValueError(f"Source experiment {expname1} does not exist.")
    if os.path.exists(target_dir):
        if clean:
            print(f"Cleaning up existing target directory {target_dir}.")
            shutil.rmtree(target_dir)
        else:
            raise ValueError(f"Target experiment {expname2} already exists.")

    # Copy the source directory to the target directory
    shutil.copytree(source_dir, target_dir)

    # Remove specific files from the target directory
    for file_to_remove in [f"{expname1}.log", "sbatch.tmp.yml"]:
        file_path = os.path.join(target_dir, file_to_remove)
        if os.path.exists(file_path):
            os.remove(file_path)

    # Rename the YAML configuration file
    source_yaml = os.path.join(target_dir, f"{expname1}.yml")
    target_yaml = os.path.join(target_dir, f"{expname2}.yml")
    if os.path.exists(source_yaml):
        os.rename(source_yaml, target_yaml)

    # Update the experiment ID in the YAML file and launch script
    with open(target_yaml, "r") as f:
        yaml_content = f.read().replace(expname1, expname2)
    with open(target_yaml, "w") as f:
        f.write(yaml_content)

    launch_script = os.path.join(target_dir, "launch.sh")
    if os.path.exists(launch_script):
        with open(launch_script, "r") as f:
            script_content = f.read().replace(expname1, expname2)
        with open(launch_script, "w") as f:
            f.write(script_content)

    # Remove temporary and backup files
    for item in os.listdir(target_dir):
        if item.startswith("ecs") or item.endswith("~"):
            shutil.rmtree(os.path.join(target_dir, item), ignore_errors=True)

    print(f"Duplicated job {expname1} to {expname2}.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Duplicate job configuration for experiments.")
    parser.add_argument("-i", "--expname1", type=str, required=True, help="Source experiment name (e.g., aa00).")
    parser.add_argument("-o", "--expname2", type=str, required=True, help="Target experiment name (e.g., bb00).")
    parser.add_argument("-c", "--config", type=str, default="config.yml", help="Path to the configuration file (default: config.yml).")
    parser.add_argument("--clean", action="store_true", help="Clean up temporary files after duplication.")

    args = parser.parse_args()
    if len(args.expname1) != 4 or len(args.expname2) != 4:
        raise ValueError("Experiment names must be 4 characters long.")

    duplicate_job(args.expname1, args.expname2, args.config, clean=args.clean)
