#!/usr/bin/env python3
"""
Simple script to modify templates of namelists & scriptlib for EC-Earth4 experiments.

Usage:
    python modify-namelists.py -e <experiment_name> -t <task_name>

    -c, --config <config_file>     Path to the configuration file (default: config.yml).
    -e, --expname <experiment_name> Name of the experiment (e.g., aa00).
    -t, --task <task_name>          Task to perform. Available tasks:
                                    - disable_wave_model: Disables the wave model in the namelist.oifs.j2 file
                                      and adds 'ignore_not_found: true' to scriptlib/post-oifs.yml.

"""

import os
import sys
import argparse
from ruamel.yaml.scalarstring import PlainScalarString
from ruamel.yaml.comments import CommentedMap
from yaml_util import load_yaml, save_yaml
from yaml_util import noparse_block, list_block

# OIFS namelist modifications

def disable_wave_model(expname, config):
    """
    Modifies the namelist.oifs.j2 file to disable the wave model 
    and add 'ignore_not_found: true' to scriptlib/post-oifs.yml.

    Args:
        namelist_path (str): Path to the namelist.oifs.j2 file.
    """

    job_dir = os.path.join(config['job_dir'], expname)
    namelist_path = os.path.join(job_dir, "templates", "oifs", "namelist.oifs.j2")
    if not os.path.exists(namelist_path):
        print(f"Error: {namelist_path} does not exist.")
        return

    try:
        with open(namelist_path, 'r') as file:
            lines = file.readlines()

        with open(namelist_path, 'w') as file:
            for line in lines:
                if "LWCOU=" in line:
                    file.write("    LWCOU=false,\n")
                elif "LWCOU2W=" in line:
                    file.write("    LWCOU2W=false,\n")
                else:
                    file.write(line)

        print(f"Successfully modified {namelist_path}.")
    except Exception as e:
        print(f"An error occurred while modifying {namelist_path}: {e}")

    # Add 'ignore_not_found: true' to scriptlib/post-oifs.yml
    if not os.path.exists(namelist_path):
        print(f"Error: {namelist_path} does not exist.")
        return

    try:
        with open(namelist_path, 'r') as file:
            lines = file.readlines()

        with open(namelist_path, 'w') as file:
            for line in lines:
                file.write(line)
                if "dst: 'restart/{{\"%03d\" % (experiment.schedule.leg.num|int+1)}}/{{wamfile_new}}'" in line:
                    file.write("        ignore_not_found: true\n")

        print(f"Successfully modified {namelist_path}.")
    except Exception as e:
        print(f"An error occurred while modifying {namelist_path}: {e}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Modify namelists for EC-Earth4 experiments.")
    parser.add_argument("-e", "--expname", required=True, help="Name of the experiment (e.g., aa00).")
    parser.add_argument("-t", "--task", required=True, choices=["disable_wave_model"],
                        help="Task to perform. Available tasks:\n"
                             "- disable_wave_model: Disables the wave model in the namelist.oifs.j2 file "
                             "and adds 'ignore_not_found: true' to scriptlib/post-oifs.yml.")

    args = parser.parse_args()
    if len(args.expname) != 4:
        raise ValueError("Experiment name must be 4 characters long.")

    # load configuration file
    config = load_yaml(args.config, expand_env=True)

    # Execute the selected task
    if args.task == "disable_wave_model":
        disable_wave_model(args.expname, config)
    else:
        print(f"Error: Unknown task '{args.task}'.")
