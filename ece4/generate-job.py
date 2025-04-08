#!/usr/bin/env python3
"""
Simple tool to generate a new EC-Earth4 experiment. This combines editing the default configuration file
with the generation of a new experiment name, creating a new folder and preparing a "launch" script.
Some caution is required to edit the YAML CommentedMap.
"""

import os
import shutil
import argparse
from ruamel.yaml.scalarstring import PlainScalarString
from yaml_util import load_yaml, save_yaml
from yaml_util import noparse_block, list_block

def create_folder(expname, config, clean=False):
    """
    Create a new folder for the experiment.

    Args:
        expname (str): Name of the experiment.
        config (str): Path to the configuration file.
        clean (bool): Whether to clean up the existing folder if it exists.
    """
    # create the folder
    conf = load_yaml(config)
    job_dir = os.path.join(conf['job_dir'], expname)

    # if clean is True, remove the existing folder
    if clean:
        if os.path.exists(job_dir):
            shutil.rmtree(job_dir)
            print(f"Removed existing job directory: {job_dir}")
        else:
            print(f"No existing job directory to remove: {job_dir}")

    # check if the folder already exists
    if os.path.exists(job_dir):
        raise ValueError(f"Experiment {expname} already exists. Please choose a different name.")
    os.makedirs(job_dir, exist_ok=True)

    # copy the template files
    base_dir = os.path.join(conf["ece_dir"], "scripts", "runtime")
    for directory in ["scriptlib", "templates"]:
        shutil.copytree(os.path.join(base_dir, directory), os.path.join(job_dir, directory), dirs_exist_ok=True)
    
    print(f"Created job directory: {job_dir}")

def create_launch(expname, config):
    """
    Create a launch bash script for the experiment.

    Args:
        expname (str): Name of the experiment.
        config (str): Path to the configuration file.
    """

    conf = load_yaml(config)
    job_dir = os.path.join(conf['job_dir'], expname)
    ece_dir = conf["ece_dir"]
    platform = conf["platform"]

    bash_script = f"""#!/bin/bash

# basic script to run the ECE4 job
platform={ece_dir}/scripts/platforms/{platform}

se user-config.yml {expname}.yml ${{platform}} scriptlib/main.yml --loglevel info
"""

    # write the file
    script_path = os.path.join(job_dir, f"launch.sh")
    with open(script_path, "w") as f:
        f.write(bash_script)

    # Make it executable (if on Unix)
    os.chmod(script_path, 0o755)

    print(f"Bash script written to: {script_path}")

def generate_user_config(expname, config):

    """
    Generate a user configuration file for the experiment.
    """
    # load configuration file 
    conf = load_yaml(config)

    # define configuration
    user_config = load_yaml(os.path.join(conf["ece_dir"], "scripts", "runtime", "user-config-example.yml"))
    user_config[0]['base.context']['experiment']['run_dir'] = noparse_block(conf['run_dir']+"/{{experiment.id}}")
    user_config[0]['base.context']['experiment']['ini_dir'] =  noparse_block(conf['ini_dir'])
    user_config[0]['base.context']['experiment']['base_dir'] =  noparse_block(conf['ece_dir'])

    # save the user configuration file
    job_dir = os.path.join(conf['job_dir'], expname)
    user_config_file = os.path.join(job_dir, "user-config.yml")
    save_yaml(user_config_file, user_config)
    print(f"User configuration file written to: {user_config_file}")


def generate_job(kind, config, expname):

    # load configuration file and setup core variables
    conf = load_yaml(config)
    job_dir = os.path.join(conf['job_dir'], expname)
    exp_base_file = os.path.join(conf["ece_dir"], "scripts", "runtime", "experiment-config-example.yml")
    account = conf["account"]

    # load base template experiment
    exp_base = load_yaml(exp_base_file)
    context = exp_base[0]['base.context']
    context['experiment']['id'] = expname

    # avoid case sensitivity
    kind = kind.upper()

    # Set the experiment name
    if kind == 'AMIP':        
        context['model_config']['components'] = list_block(['oifs', 'amipfr', 'xios', 'oasis'])
        context['model_config']['oifs']['grid'] = noparse_block("{{model_config.oifs.all_grids."+conf["resolution"]["oifs"]+"}}")
        del context['model_config']['nemo']
    elif kind == 'CPLD':
        context['model_config']['components'] = list_block(['oifs', 'nemo', 'rnfm', 'xios', 'oasis'])
        context['model_config']['oifs']['grid'] = noparse_block("{{model_config.oifs.all_grids."+conf["resolution"]["oifs"]+"}}")
        context['model_config']['nemo']['grid'] = noparse_block("{{model_config.nemo.all_grids."+conf["resolution"]["nemo"]+"}}")
    
    # setup job block
    context['job']['launch']['method'] = PlainScalarString('slurm-wrapper-taskset')
    context['job']['slurm']['sbatch']['opts']['account']= account
    context['job']['slurm']['sbatch']['opts']['qos'] = 'np'
    context['job']['slurm']['sbatch']['opts']['time'] = 180
    context['job']['slurm']['sbatch']['opts']['ntasks-per-core'] = 1

    # delete the not wrapper tasket block: this might change in the future
    del exp_base[1]

    # default one node configuration for AMIP
    if kind == "AMIP":
        exp_base[1]['base.context']['job']['groups'] = [{'nodes': 1, 'xios': 1, 'oifs': 125, 'amipfr': 1}]
    # default two node configuration for CPLD
    elif kind == "CPLD":
        exp_base[1]['base.context']['job']['groups'] = [
            { 'nodes': 1, 'xios': 1, 'oifs': 126, 'rnfm': 1 },
            { 'nodes': 1, 'oifs': 49, 'nemo': 77 },
        ]

    # write the file
    yaml_path = os.path.join(job_dir, f'{expname}.yml')
    save_yaml(os.path.join(job_dir, f'{expname}.yml'), exp_base)
    print(f"YAML script script written to: {yaml_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate job configuration for experiments.")
    parser.add_argument("-k", "--kind", type=str, required=True, help="Type of experiment (e.g., AMIP).")
    parser.add_argument("-c","--config", type=str, help="Type of experiment (e.g., AMIP).", default="config.yml")
    parser.add_argument("-e", "--expname", type=str, required=True, help="Experiment name (e.g., aa00).")
    parser.add_argument("--clean", action="store_true", help="Clean up the experiment folder.")

    args = parser.parse_args()
    if args.kind.upper() not in ["AMIP", "CPLD"]:
        raise ValueError("Invalid experiment type. Choose either 'AMIP' or 'CPLD'.")
    if len(args.expname) != 4:
        raise ValueError("Experiment name must be 4 characters long.")
    create_folder(args.expname, args.config, args.clean)
    generate_job(args.kind, args.config, args.expname)
    generate_user_config(args.expname, args.config)
    create_launch(args.expname, args.config)