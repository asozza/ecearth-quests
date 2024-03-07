"""
Module to load and manipulate yaml file with ruamel package

Authors
Matteo Nurisso (CNR-ISAC, Mar 2024)
"""
import os
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedSeq

def load_yaml(file: str= None, ruamel_type: str = 'rt'):
    """
    Load yaml file with ruamel.yaml package

    Args:
        file (str): a file path to the yaml
        ruamel_type (str, optional): the type of YAML initialisation.
                                     Default is 'rt' (round-trip)

    Returns:
        A dictionary with the yaml file keys
    """

    if not os.path.exists(file):
        raise ValueError(f'File {file} not found: you need to have this configuration file!')

    yaml = YAML(typ=ruamel_type)

    # Load the YAML file as a text string
    with open(file, 'r', encoding='utf-8') as file:
        yaml_text = file.read()
    
    cfg = yaml.load(yaml_text)

    if isinstance(cfg, CommentedSeq):
        cfg = cfg[0]

    return cfg