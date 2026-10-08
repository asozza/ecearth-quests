#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a command line tool to create with CDO the descriptor file necessary to
perform a wide range of horizontal interpolation of initial conditions

Authors
Paolo Davini (CNR-ISAC, Nov 2023)
"""

import os
import subprocess
from utils import extract_grid_info, ecmwf_grid

# grid list
grids = ['TL63', 'TL95', 'TCO95', 'TL159', 'TCO199', 'TCO319', 'TCO399']

# this runs only on Atos if you have the rights to read this folder
OIFS_BC="/lus/h2resw01/fws1/mirror/lb/project/rdxdata/climate/climate.v015"


for grid in grids:

    grid = grid.upper()
    print('Processing ' + grid + '...')

    grid_type, truncation, _ = extract_grid_info(grid + 'L31')
    ecmwf_kind = str(truncation) + ecmwf_grid(grid_type)
    REF_FILE='10_bats_glcc.grb'
    file_path = os.path.join(OIFS_BC, ecmwf_kind, REF_FILE)
    target_path = os.path.join('grids', grid + '.txt')
    #pippo = cdo.griddes(input=file_path, stdout=target_path)
    with open(target_path, 'w', encoding='utf8') as fff:
        subprocess.run(['cdo', 'griddes', file_path], stdout=fff, check=True)

