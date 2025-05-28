#!/usr/bin/env python3

import os
import shutil
import argparse

from yaml_util import load_yaml, save_yaml

import generate_job as gj

########################################################################

def duplicate_and_perturb(parnam_ok, parval_ok, expname1, expname2, config, clean=False):
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

    # Change context with new tuning param
    exp_base = load_yaml(target_yaml)
    context = exp_base[0]['base.context']
    context['experiment']['id'] = expname2

    namelists = ['namcumf', 'namcldp']
    parlist = dict()
    parlist['namcumf'] = ['RPRCON', 'ENTRORG', 'DETRPEN', 'ENTRDD', 'RMFDEPS']
    parlist['namcldp'] = 'RVICE RLCRITSNOW RSNOWLIN2 RCLDIFF RCLDIFF_CONVI'.split()
 
    for namls in namelists:
        if parnam_ok in parlist[namls]:
            print(f'Changing {parnam_ok} to: {parval_ok}')
            context['model_config']['oifs']['tuning'][namls][parnam_ok] = parval_ok

    # write the file
    yaml_path = os.path.join(target_dir, f'{expname2}.yml')
    save_yaml(os.path.join(target_dir, f'{expname2}.yml'), exp_base)
    print(f"YAML script script written to: {yaml_path}")

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

    print(f"Duplicated job {expname1} to {expname2} with a change of {parnam_ok}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Duplicate job configuration for experiments.")
    parser.add_argument("-i", "--expname1", type=str, required=True, help="Source experiment name (e.g., aa00).")
    parser.add_argument("-c", "--config", type=str, default="config.yml", help="Path to the configuration file (default: config.yml).")
    parser.add_argument("--clean", action="store_true", help="Clean up temporary files after duplication.")

    args = parser.parse_args()
    if len(args.expname1) != 4:
        raise ValueError("Experiment names must be 4 characters long.")

    params = load_yaml("tuning_params.yml")

    basename = args.expname1[:2]
    icount = 1

    for par in params:
        # perturb in both direction
        for pert in [-20, 20]:
            new_val = round(params[par] * (1 + pert/100.), 10)
            expnam = basename + f"{icount:02d}"

            duplicate_and_perturb(par, new_val, args.expname1, expnam, args.config, clean=args.clean)

            icount +=1