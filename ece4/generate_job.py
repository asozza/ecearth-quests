#!/usr/bin/env python3
"""
Simple tool to generate a new EC-Earth4 experiment. This combines editing the default configuration file
with the generation of a new experiment name, creating a new folder and preparing a "launch" script.
Some caution is required to edit the YAML CommentedMap and CommentedSeq objects.

Usage:
    python generate-job.py -k <experiment_type> -e <experiment_name> [-c <config_file>] [--clean]

    -k, --kind <experiment_type>   Type of experiment (e.g., AMIP, CPLD, OMIP).
    -e, --expname <experiment_name> Name of the experiment (e.g., aa00).
    -c, --config <config_file>     Path to the configuration file (default: config.yml).
    --clean                        Clean up the experiment folder if it exists.
"""

import os
import shutil
import argparse
from ruamel.yaml.scalarstring import PlainScalarString
from ruamel.yaml.comments import CommentedMap
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
   
    job_dir = os.path.join(config['job_dir'], expname)

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
    base_dir = os.path.join(config["ece_dir"], "scripts", "runtime")
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

    job_dir = os.path.join(config['job_dir'], expname)
    ece_dir = config["ece_dir"]
    platform = config["platform"]

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

    Args:
        kind (str): Type of experiment (e.g., AMIP).
        expname (str): Name of the experiment.
        config (str): Path to the configuration file.
    """

    # load configuration file 
    src_dir = config['ece_dir']

    # define configuration
    user_config = load_yaml(os.path.join(src_dir, "scripts", "runtime", "user-config-example.yml"))
    user_config[0]['base.context']['experiment']['run_dir'] = noparse_block(config['run_dir']+"/{{experiment.id}}")
    user_config[0]['base.context']['experiment']['ini_dir'] =  noparse_block(config['ini_dir'])
    user_config[0]['base.context']['experiment']['base_dir'] =  noparse_block(src_dir)

    # save the user configuration file
    job_dir = os.path.join(config['job_dir'], expname)
    user_config_file = os.path.join(job_dir, "user-config.yml")
    save_yaml(user_config_file, user_config)
    print(f"User configuration file written to: {user_config_file}")


def add_tuning_section(config, context):
    """
    Adds custom tuning parameters directly into the jobscript.
    """

    print('Changing tuning parameters!')

    if 'tuning' not in config:
        raise ValueError('Tuning section not found in config!')
    
    if 'tuning' not in context['model_config']['oifs']:
        context['model_config']['oifs']['tuning'] = dict()

    namelists = ['namcumf', 'namcldp']
    parlist = dict()
    parlist['namcumf'] = ['RPRCON', 'ENTRORG', 'DETRPEN', 'ENTRDD', 'RMFDEPS']
    parlist['namcldp'] = 'RVICE RLCRITSNOW RSNOWLIN2 RCLDIFF RCLDIFF_CONVI'.split()
 
    params = config['tuning']

    for namls in namelists:
        if any([par in params for par in parlist[namls]]):
            context['model_config']['oifs']['tuning'][namls] = dict()

        for parnam in parlist[namls]:
            if parnam in params:
                print(f'Changing {parnam} to: {params[parnam]}')
                context['model_config']['oifs']['tuning'][namls][parnam] = params[parnam]
    
    return context


def generate_job(kind, config, expname):
    """
    Generate a job configuration file for the experiment.

    Args:
        kind (str): Type of experiment (e.g., AMIP).
        config (str): Path to the configuration file.
        expname (str): Name of the experiment.
    """

    # load configuration file and setup core variables
    job_dir = os.path.join(config['job_dir'], expname)
    src_dir = config['ece_dir']
    exp_base_file = os.path.join(src_dir, "scripts", "runtime", "experiment-config-example.yml")

    # load base template experiment
    exp_base = load_yaml(exp_base_file)
    context = exp_base[0]['base.context']
    context['experiment']['id'] = expname

    # avoid case sensitivity
    kind = kind.upper()

    # adjust schedule settings
    startdate = config['schedule']['start']
    enddate = config['schedule']['end']
    freq = config['schedule']['freq']
    context['experiment']['schedule']['all'].value = f"DTSTART:{startdate} RRULE:FREQ={freq};UNTIL={enddate}"
    context['experiment']['schedule']['nlegs'] = config['schedule']['nlegs']

    # Set the experiment name
    if kind == 'AMIP':        
        context['model_config']['components'] = list_block(['oifs', 'amipfr', 'xios', 'oasis'])
        context['model_config']['oifs']['grid'] = noparse_block("{{model_config.oifs.all_grids."+config["resolution"]["oifs"]+"}}")
        del context['model_config']['nemo']
    elif kind == 'CPLD':
        context['model_config']['components'] = list_block(['oifs', 'nemo', 'rnfm', 'xios', 'oasis'])
        context['model_config']['oifs']['grid'] = noparse_block("{{model_config.oifs.all_grids."+config["resolution"]["oifs"]+"}}")
        context['model_config']['nemo']['grid'] = noparse_block("{{model_config.nemo.all_grids."+config["resolution"]["nemo"]+"}}")
    elif kind == 'OMIP':
        context['model_config']['components'] = list_block(['nemo', 'xios'])
        context['model_config']['nemo']['grid'] = noparse_block("{{model_config.nemo.all_grids."+config["resolution"]["nemo"]+"}}")
        del context['model_config']['oifs']
        del context['model_config']['oasis']

    if 'tuning' in config:
        context = add_tuning_section(config, context)
    
    if 'forcing' in config:
        if 'cmip' not in context['experiment']['forcing'] or context['experiment']['forcing']['cmip'] is None:
            context['experiment']['forcing']['cmip'] = dict()
            context['experiment']['forcing']['cmip']['fixyear'] = -1

        if config['forcing'] == 'historical':
            context['experiment']['forcing']['cmip']['fixyear'] = -1
        elif config['forcing'] == 'preindustrial':
            context['experiment']['forcing']['cmip']['fixyear'] = 1850
        elif isinstance(config['forcing'], int):
            context['experiment']['forcing']['cmip']['fixyear'] = config['forcing']
        else:
            raise ValueError('forcing not understood, use one of preindustrial, historical or a custom year')

    # setup job block
    context['job']['launch']['method'] = PlainScalarString(config['launch-method'])
    if config['launch-method'] != 'slurm-wrapper-taskset':
        context['job']['launch']['shell'] = CommentedMap()
        context['job']['launch']['shell']['script'] = PlainScalarString('run-srun-multiprog.sh')

    context['job']['slurm']['sbatch']['opts']['account']= config["account"]
    context['job']['slurm']['sbatch']['opts']['qos'] = 'np'
    context['job']['slurm']['sbatch']['opts']['time'] = 180
    context['job']['slurm']['sbatch']['opts']['ntasks-per-core'] = 1
    
    if config['launch-method'] == 'slurm-wrapper-taskset':
        # delete the not wrapper-tasket block (this might change in the future)
        del exp_base[1]
        exp_base.yaml_set_comment_before_after_key(1, before='\n')

        # default one node configuration for AMIP
        if kind == "AMIP":
            exp_base[1]['base.context']['job']['groups'] = [{'nodes': 1, 'xios': 1, 'oifs': 126, 'amipfr': 1}]
        # default two node configuration for CPLD
        elif kind == "CPLD":
            exp_base[1]['base.context']['job']['groups'] = [
                { 'nodes': 1, 'xios': 1, 'oifs': 126, 'rnfm': 1 },
                { 'nodes': 1, 'oifs': 38, 'nemo': 90 },
            ]
        # default one node configuration for OMIP
        elif kind == "OMIP":
            exp_base[1]['base.context']['job']['groups'] = [
                { 'nodes': 1, 'xios': 1, 'nemo': 90 },
            ]
    else:
        # delete the wrapper-taskset block
        del exp_base[2]

        # default configurations for AMIP and CPLD
        if kind == "AMIP":
            exp_base[1]['base.context']['job'] = {
                'oifs': {'ntasks': 126, 'ntasks_per_node': 128},
                'xios': {'ntasks': 1, 'ntasks_per_node': 1},
                'amipfr': {'ntasks': 1, 'ntasks_per_node': 1}
            }
        elif kind == "CPLD":
            exp_base[1]['base.context']['job'] = {
                'oifs': {'ntasks': 164, 'ntasks_per_node': 128},
                'nemo': {'ntasks': 90, 'ntasks_per_node': 128},            
                'xios': {'ntasks': 1, 'ntasks_per_node': 1},
                'rnfm': {'ntasks': 1, 'ntasks_per_node': 1}
            }
        elif kind == "OMIP":
            exp_base[1]['base.context']['job'] = {
                'nemo': {'ntasks': 90, 'ntasks_per_node': 128},
                'xios': {'ntasks': 1, 'ntasks_per_node': 1}
            }

    # write the file
    yaml_path = os.path.join(job_dir, f'{expname}.yml')
    save_yaml(os.path.join(job_dir, f'{expname}.yml'), exp_base)
    print(f"YAML script script written to: {yaml_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate job configuration for experiments.")
    parser.add_argument("-k", "--kind", type=str, required=True, help="Type of experiment (e.g., AMIP, CPLD, OMIP).")
    parser.add_argument("-c","--config", type=str, help="YAML configuration file", default="config.yml")
    parser.add_argument("-e", "--expname", type=str, required=True, help="Experiment name (e.g., aa00).")
    parser.add_argument("--clean", action="store_true", help="Clean up the experiment folder.")

    args = parser.parse_args()
    if args.kind.upper() not in ["AMIP", "CPLD", "OMIP"]:
        raise ValueError("Invalid experiment type. Choose either 'AMIP', 'CPLD' or 'OMIP'.")
    if len(args.expname) != 4:
        raise ValueError("Experiment name must be 4 characters long.")
    
    # load configuration file
    config = load_yaml(args.config, expand_env=True)

    create_folder(args.expname, config, args.clean)
    generate_job(args.kind, config, args.expname)
    generate_user_config(args.expname, config)
    create_launch(args.expname, config)
