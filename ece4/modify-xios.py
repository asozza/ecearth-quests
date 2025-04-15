#!/usr/bin/env python3
"""
Simple script to modify XIOS XML files for EC-Earth4 experiments.
This script adds 'append="true"' to every <file> tag in the given XML file and removes <file> tags 
with name_suffix containing a specific keyword.

Usage:
    ./modify-xios.py <xml_file> <keywords>
    
    - xml_file: Path to the XML file to be modified.
    - keywords: List of keywords to search for in the name_suffix attribute of <file> tags to be removed.
"""


import os
import shutil
import argparse
import xml.etree.ElementTree as ET

def add_append_to_file_tags(xml_file, backup=True):
    """ 
    Adds 'append="true"' to every <file> tag in the given XML file.
    """

    # Backup
    if backup:
        backup_file = xml_file + ".bak"
        if not os.path.exists(backup_file):
            os.system(f"cp {xml_file} {backup_file}")

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Add 'append="true"' to every <file> tag
    for file_tag in root.iter('file'):
        if 'append' not in file_tag.attrib:
            file_tag.set('append', 'true')

    tree.write(xml_file, encoding='utf-8', xml_declaration=True)
    print(f"Updated: {xml_file}")


def remove_file_by_suffix(xml_file, keyword, backup=True):
    """
    Removes <file> tags with name_suffix containing the keyword from the given XML file.
    """

    # Backup
    if backup:
        backup_file = xml_file + ".bak"
        if not os.path.exists(backup_file):
            os.system(f"cp {xml_file} {backup_file}")

    tree = ET.parse(xml_file)
    root = tree.getroot()

    files_to_remove = []

    # Trova tutti i blocchi <file> con name_suffix contenente la keyword
    for parent in root.iter():
        for file_tag in list(parent):
            if file_tag.tag == 'file':
                name_suffix = file_tag.attrib.get('name_suffix', '')
                if keyword in name_suffix:
                    files_to_remove.append((parent, file_tag))

    # Rimuove i blocchi trovati
    for parent, file_tag in files_to_remove:
        parent.remove(file_tag)
        print(f"Removed file tag with name_suffix='{file_tag.attrib.get('name_suffix')}'")

    tree.write(xml_file, encoding='utf-8', xml_declaration=True)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Modify XIOS XML files for EC-Earth4 experiments.")
    parser.add_argument("xml_file", type=str, help="Path to the XML file to be modified.")
    parser.add_argument("keywords", nargs='+', help="List of keywords for name_suffix tags.")

    args = parser.parse_args()
    # Add 'append="true"' to <file> tags
    add_append_to_file_tags(args.xml_file)

    # Remove <file> tags with specified keywords (e.g. '_1m_diaptr2d', '_1m_diaptr3d')
    for keyword in args.keywords:
        remove_file_by_suffix(args.xml_file, keyword)
        
