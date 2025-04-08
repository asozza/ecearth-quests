"""
Module to load and manipulate yaml file with ruamel package

Authors
Matteo Nurisso (CNR-ISAC, Mar 2024)
"""
import os
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import PlainScalarString
from ruamel.yaml.comments import TaggedScalar, CommentedSeq

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

    return cfg


def save_yaml(path: str = None, cfg: dict = None, ruamel_type: str = 'rt'):
    """
    Save dictionary to a yaml file with ruamel.yaml package

    Args:
        path (str): a file path to the yaml
        cfg (dict): a dictionary to be dumped
        ruamel_type (str, optional): the type of YAML initialisation.
                                    Default is 'rt' (round-trip)
    """
    # Initialize YAML object
    yaml = YAML(typ=ruamel_type)

    # Check input
    if path is None:
        raise ValueError('File not defined')
    if cfg is None:
        raise ValueError('Content cfg not defined')

    # Dump to file
    with open(path, 'w', encoding='utf-8') as path:
        yaml.dump(cfg, path)

    return None

def noparse_block(value):
    """
    Create a block scalar with the !noparse tag.
    """
    return TaggedScalar(value, tag="!noparse", style='"')

def list_block(value):
    """
    Create a PlanScalar with a list
    """
    return list_compact([PlainScalarString(x) for x in value])

def list_compact(list):

    """
    Create a compact list in CommentedSeq format.
    """
    comm = CommentedSeq(list)
    comm.fa.set_flow_style()
    return comm

