#!/usr/bin/env python3
"""
Simple script to perturbate a job configuration for tuning experiments.
This script duplicates an existing job configuration, modifies a specified parameter, and prepares the new job for execution.
This is useful for generating multiple tuning experiments based on a reference configuration.
Usage:
    python create_tuning_ensemble.py -e <experiment_name> [-c <config_file>] [--clean]

    -e, --expname <experiment_name> Name of the reference experiment (e.g., aa00).
    -c, --config <config_file>     Path to the configuration file (default: config.yml).
    --clean                        Clean up the experiment folder if it exists.    
""" 

import os
import shutil
import argparse

from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import PlainScalarString
from ruamel.yaml.comments import CommentedMap
from yaml_util import load_yaml, save_yaml
from yaml_util import noparse_block, list_block

yaml = YAML()
yaml.preserve_quotes = True

#########################################################################
# Global variables

components = ['oifs', 'nemo']

########################################################################

def perturbate(parname, parvalue, expname1, expname2, config, clean=False, skip_existing = False):
    """
    Duplicate an existing job to create a new one with a different experiment ID 
    and perturb a specified parameter.

    Args:
        expname1 (str): Name of the source experiment.
        expname2 (str): Name of the target experiment.
        config (str): Path to the configuration file.
    """

    # load the configuration file
    job_dir = config['job_dir']
    source_dir = os.path.join(job_dir, expname1)
    target_dir = os.path.join(job_dir, expname2)

    if not os.path.exists(source_dir):
        raise ValueError(f"Reference experiment {expname1} does not exist.")
    if os.path.exists(target_dir):
        if clean:
            print(f"Cleaning up existing target directory {target_dir}.")
            shutil.rmtree(target_dir)
        else:
            if skip_existing:
                print(f'Job already exists!')
                print(f"Leaving {expname2} with a change of {parname} as is, going on...")
                return
            else:
                raise ValueError(f"Target experiment {expname2} already exists.")

    # Copy the source directory to the target directory
    shutil.copytree(source_dir, target_dir)

    # Remove specific files from the target directory
    for file_to_remove in [f"{expname1}.log", "sbatch.tmp.yml"]:
        file_path = os.path.join(target_dir, file_to_remove)
        if os.path.exists(file_path):
            os.remove(file_path)

    # Rename YAML configuration file
    source_yaml = os.path.join(target_dir, f"{expname1}.yml")
    target_yaml = os.path.join(target_dir, f"{expname2}.yml")
    if os.path.exists(source_yaml):
        os.rename(source_yaml, target_yaml)

    # Change YAML configuration file
    exp_base = load_yaml(target_yaml)
    context = exp_base[0]['base.context']
    context['experiment']['id'] = expname2
    context['model_config']['tuning_file'] = noparse_block("{{se.cli.cwd}}/templates/tuning-example.yml") # enable tuning file
    save_yaml(target_yaml, exp_base) # write modified YAML file
    print(f"YAML configuration file written to: {target_yaml}")

    ################################################################

    # Update tuning file
    tuning_yaml = os.path.join(target_dir, "templates", "tuning-example.yml")
    tuning_base = load_yaml(tuning_yaml)
    tuning_context = tuning_base[0]['base.context']

    print(f"Updating {comp}.{namelist}.{parname} to {parvalue}")
    if comp not in tuning_context['model_config']:
        tuning_context['model_config'][comp] = dict()
        tuning_context['model_config'][comp]['tuning'] = dict()
        tuning_context['model_config'][comp]['tuning'][namelist] = dict()
    tuning_context['model_config'][comp]['tuning'][namelist][parname] = parvalue
    save_yaml(tuning_yaml, tuning_base)
    print(f"tuning YAML file written to: {tuning_yaml}")

    # Update the launch script
    lfile = os.path.join(target_dir, "launch.sh")
    if os.path.exists(lfile):
        with open(lfile, "r") as f:
            content = f.read().replace(expname1, expname2)
        with open(lfile, "w") as f:
            f.write(content)

    # Remove temporary and backup files
    for item in os.listdir(target_dir):
        if item.startswith("ecs") or item.endswith("~"):
            shutil.rmtree(os.path.join(target_dir, item), ignore_errors=True)

    print(f"Duplicated job {expname1} to {expname2} with a change of {parname}.")

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Perturbate job configurations for tuning experiments.")
    parser.add_argument("-e", "--expname", type=str, required=True, help="Reference experiment name (e.g., aa00).")
    parser.add_argument("-p", "--percent", type=float, default=20.0, help="Percentage to perturb the parameter (default: 20.0).")
    parser.add_argument("-c", "--config", type=str, default="config.yml", help="Path to the configuration file (default: config.yml).")
    parser.add_argument("-t", "--tuning_params", type=str, default="tuning_params.yml", help="Path to the tuning file with all params to perturb (default: tuning_params.yml).")
    parser.add_argument("--clean", action="store_true", help="Clean up temporary files after duplication.")
    parser.add_argument("--skip_existing", action="store_true", help="Leave existing exps where they are (if adding new params)")

    args = parser.parse_args()
    if len(args.expname) != 4:
        raise ValueError("Experiment names must be 4 characters long.")

    params = load_yaml(args.tuning_params)
    config = load_yaml(args.config, expand_env=True)

    percent = args.percent
    icount = 1

    # cycle through all components and namelists
    for comp in components:
        if comp in params:
            print(f"Perturbating component: {comp}")
        else:
            print(f"Component {comp} not found in tuning parameters, skipping.")
            continue        
        for namelist in params[comp]:
            for parname, parvalue in params[comp][namelist].items():

                # Apply perturbations in both directions
                for perc in [-percent, percent]:
                    new_value = round(float(parvalue) * (1 + perc / 100.0), 10)
                    expname2 = args.expname[:2] + f"{icount:02d}"
                    perturbate(parname, new_value, args.expname, expname2, config, clean=args.clean, skip_existing=args.skip_existing)
                    icount += 1

