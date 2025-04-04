#!/usr/bin/env python3

import argparse
from yaml_util import load_yaml, save_yaml
from ruamel.yaml.scalarstring import PlainScalarString
from ruamel.yaml.comments import TaggedScalar
import os

def generate_job(kind, config, expname):

    conf = load_yaml(config)

    exp_base_file = os.path.join(conf["ece_dir"], "scripts", "runtime", "experiment-config-example.yml")
    exp_base = load_yaml(exp_base_file)
    context = exp_base[0]['base.context']

    # Set the experiment name
    if kind == 'AMIP':
        context['experiment']['id'] = expname
        context['model_config']['components'] = [PlainScalarString(x) for x in ['oifs', 'amipfr', 'xios', 'oasis']]
        val = TaggedScalar("{{model_config.oifs.all_grids."+conf["resolution"]["oifs"]+"}}", tag="!noparse")
        context['model_config']['oifs']['grid'] = val
        del context['model_config']['nemo']
    
    context['job']['launch']['method'] = PlainScalarString('slurm-wrapper-taskset')
    del exp_base[1]

    if kind == "AMIP":
        exp_base[1]['base.context']['job']['groups'] = [{'nodes': 1, 'xios': 1, 'oifs': 125, 'amipfr': 1}]

    save_yaml(f'{expname}.yml', exp_base)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate job configuration for experiments.")
    parser.add_argument("--kind", type=str, required=True, help="Type of experiment (e.g., AMIP).")
    parser.add_argument("--config", type=str, help="Type of experiment (e.g., AMIP).", default="config.yml")
    parser.add_argument("--expname", type=str, required=True, help="Experiment name (e.g., aa00).")

    args = parser.parse_args()
    generate_job(args.kind, args.config, args.expname)