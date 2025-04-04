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
from ruamel.yaml.comments import TaggedScalar
from yaml_util import load_yaml, save_yaml

def create_folder(expname, config):
    """
    Create a new folder for the experiment.
    """
    # create the folder
    conf = load_yaml(config)
    exp_dir = os.path.join(conf['job_dir'], expname)
    if os.path.exists(exp_dir):
        raise ValueError(f"Experiment {expname} already exists. Please choose a different name.")
    os.makedirs(exp_dir, exist_ok=True)

    base_dir = os.path.join(conf["ece_dir"], "scripts", "runtime")
    for directory in ["scriptlib", "templates"]:
        shutil.copytree(os.path.join(base_dir, directory), os.path.join(exp_dir, directory), dirs_exist_ok=True)


def generate_job(kind, config, expname):

    # load configuration file and setup core variables
    conf = load_yaml(config)
    exp_dir = os.path.join(conf['job_dir'], expname)
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

    # delete the not wrapper tasket block
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

    save_yaml(os.path.join(exp_dir, f'{expname}.yml'), exp_base)

def noparse_block(value):
    """
    Create a block scalar with the !noparse tag.
    """
    return TaggedScalar(value, tag="!noparse")

def list_block(value):
    """
    Create a PlanScalar with a list
    """
    return [PlainScalarString(x) for x in value]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate job configuration for experiments.")
    parser.add_argument("--kind", type=str, required=True, help="Type of experiment (e.g., AMIP).")
    parser.add_argument("--config", type=str, help="Type of experiment (e.g., AMIP).", default="config.yml")
    parser.add_argument("--expname", type=str, required=True, help="Experiment name (e.g., aa00).")

    args = parser.parse_args()
    if args.kind.upper() not in ["AMIP", "CPLD"]:
        raise ValueError("Invalid experiment type. Choose either 'AMIP' or 'CPLD'.")
    if len(args.expname) != 4:
        raise ValueError("Experiment name must be 4 characters long.")
    create_folder(args.expname, args.config)
    generate_job(args.kind, args.config, args.expname)