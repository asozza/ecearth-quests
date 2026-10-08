#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HPC tools

Authors
Alessandro Sozza (CNR-ISAC, Mar 2024)
"""

import os
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt

########################################################
# Searching by patterns

# search for values preceded by label
def read_value_preceded_by_label(file_path, label):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                pattern = fr'{label}\s*[:=]?\s*(\d+(?:\.\d+)?)'
                match = re.search(pattern, line, re.IGNORECASE)
                if match:
                    value= float(match.group(1))
                    return value
        # if label is not found
        return None
    except FileNotFoundError:
        print(file_path," not found")
        return None

# search for values followed by label
def read_value_followed_by_label(file_path, label):
    try:
        with open(file_path, 'r') as file:
            for	line in	file:
                pattern = fr'(\d+(?:\.\d+)?)\s*[^0-9a-zA-Z]*{label}[^0-9a-zA-Z]*\s*[:=]?'
                match = re.search(pattern, line, re.IGNORECASE)
                if match:
                    value= float(match.group(1))
                    return value
        # if label is not found
        return None
    except FileNotFoundError:
        print(file_path," not found")
        return None

# search for value constrained by two labels
def read_value_constrained_by_two_labels(file_path, label):
    try:
        start_search = False
        with open(file_path, 'r') as file:
            for line in file:
                pattern = fr"Load balance analysis"
                match = re.search(pattern, line, re.IGNORECASE)
                if match:
                    start_search = True                    
                if  start_search:
                    parts = line.strip().split('/')
                    if len(parts) == 3 and parts[0].strip() == label:
                        value1 = float(parts[1].strip())
                        value2 = float(parts[2].strip())
                        return value1, value2
        # if label is not found
        return None, None
    except FileNotFoundError:
        print(file_path," not found")
        return None, None

# search for value below a label
def read_value_below_a_label(file_path, label):
    try:
        start_search = False
        with open(file_path, 'r') as file:
            for line in file:
                if start_search and line.strip().startswith(':'):
                    value = float(line.split(':')[1].strip())
                    return value
                if label in line:
                    start_search = True                    
        # if label is not found
        return None
    except FileNotFoundError:
        print(file_path," not found")
        return None

# read time (in format hh:mm:ss) preceded by label
def read_time_preceded_by_label(file_path, label):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                pattern = fr'{label}:\s*(\d+):(\d+):(\d+\.\d+)'
                match = re.search(pattern, line, re.IGNORECASE)
                if match:
                    hours = int(match.group(1))
                    minutes = int(match.group(2))
                    seconds = float(match.group(3))
                    return hours, minutes, seconds
        # if label is not found
        return None, None, None
    except FileNotFoundError:
        print(file_path, " not found")
        return None, None, None

#############################################################
# i/o operations, folders and legs
    
def mainpath():

    # read from yaml file
    with open('paths.yaml') as jf:
        config = yaml.load(jf, Loader=yaml.FullLoader)
    if 'folders' in config:
        folders = config['folders']
    for ifo in folders:
        path = ifo.get('path', '')

    return path

def folders(expname):

    path = mainpath()
    dirs = {
        'exp': os.path.join(path, expname),
        'log': os.path.join(path, expname, "log"),
    }

    return dirs

# extract final leg
def endleg(expname):

    dirs = folders(expname)
    legfile = os.path.join(dirs['exp'], 'leginfo.yml')
    with open(legfile, 'r', encoding='utf-8') as file:
        leginfo = yaml.load(file, Loader=yaml.FullLoader)
    legs = leginfo['base.context']['experiment']['schedule']['leg']['num']

    return legs

def selectleg(expname, leg):

    if leg == 0:
        legs = endleg(expname)
    else:
        legs = leg

    return legs

#############################################################
# actions: compute, group_by, collect

def compute_nodes(expname):

    init=1
    dirs = folders(expname)
    path = os.path.join(dirs['log'], str(init).zfill(3), 'NODE.001_01')
    npa = read_value_preceded_by_label(path, 'NPROC') # nprocs of oifs, npa    
    path = os.path.join(dirs['log'], str(init).zfill(3), 'ocean.output')
    npo = read_value_preceded_by_label(path, 'jpnij') # nprocs of nemo, npo

    return npa,npo

def compute_sypd(expname, legs):

    dirs = folders(expname)    
    sypd = np.array([])
    for leg in range(1,int(legs)):
        path = os.path.join(dirs['log'], str(leg).zfill(3), 'timing.log')
        value = read_value_followed_by_label(path, 'SYPD') # simulated year per day, sypd
        sypd = np.append(sypd, value)
    
    return sypd

def compute_syph(expname, legs):

    dirs = folders(expname)
    syph = np.array([])
    for leg in range(1,int(legs)):
        path = os.path.join(dirs['log'], str(leg).zfill(3), 'timing.log')
        value = read_value_followed_by_label(path, 'SYPD') # simulated year per day, sypd
        syph = np.append(syph, 24./value)
    
    return syph

def compute_chpsy(expname, legs):

    dirs = folders(expname)
    npa,npo = compute_nodes(expname)
    nptot=npa+npo+1
    chpsy = np.array([])
    for leg in range(1,int(legs)):
        path = os.path.join(dirs['log'], str(leg).zfill(3), 'timing.log')
        value = read_value_followed_by_label(path, 'SYPD') # simulated year per day, sypd
        chpsy = np.append(chpsy, 24.*nptot/value)
    
    return chpsy

def compute_elapsedtime(expname, legs):

    dirs = folders(expname)
    elapsedtime = np.array([])
    for leg in range(1,int(legs)):
        path = os.path.join(dirs['log'], str(leg).zfill(3), 'timing.log')
        hours, minutes, seconds = read_time_preceded_by_label(path, 'elapsed time')    
        value = hours*3600.+minutes*60.+seconds
        elapsedtime = np.append(elapsedtime, value)

    return elapsedtime

def compute_sbu(expname, legs):

    dirs = folders(expname)
    npa,npo = compute_nodes(expname)
    nptot=npa+npo+1
    sbu = np.array([])
    for leg in range(1,int(legs)):
        path = os.path.join(dirs['log'], str(leg).zfill(3), 'timing.log')
        hours, minutes, seconds = read_time_preceded_by_label(path, 'elapsed time')    
        time = hours*3600.+minutes*60.+seconds
        value = nptot*time*470410408./(86400.*8098.*128.) # system billing units, sbu
        sbu = np.append(sbu, value)

    return sbu

#############################################################
# collect variables from multiple runs
def group_by_nptot_sypd(expnames, leg):
    
    npo = np.array([])
    npa = np.array([])
    for expname in expnames:
        np1,np2 = compute_nodes(expname)
        npa = np.append(npa, np1)
        npo = np.append(npo, np2)    
    nptot = npo+npa+1

    ave=np.array([])
    for expname in expnames:
        legs = selectleg(expname, leg)
        sypd = compute_sypd(expname, legs)
        ave = np.append(ave,np.mean(sypd))

    return nptot,ave

def collect_sypd(expnames, leg):

    ave=np.array([])
    for expname in expnames:        
        legs = selectleg(expname, leg)
        sypd = compute_sypd(expname, legs)
        ave = np.append(ave,np.mean(sypd))

    return ave

#############################################################
# plots
# ISSUE1 --> generalize functions for any variable
# ISSUE2 --> cumulative SBU

def plot_sypd_vs_time(expname, leg):

    legs = selectleg(expname, leg)
    sypd = compute_sypd(expname, legs)
    x = [leg for leg in range(1,legs)]
    plt.xlabel("leg")
    plt.ylabel("SYPD")
    pp = plt.plot(x,sypd)

    return pp

# multiplot
# ISSUE: it can be generalized using 4 labels
def multiplot_vs_time(expname, leg):

    legs = selectleg(expname, leg)
    x = [leg for leg in range(1,legs)]
    sypd = compute_sypd(expname, legs)
    elapsedtime = compute_elapsedtime(expname, legs)
    chpsy = compute_chpsy(expname, legs)
    sbu = compute_sbu(expname, legs)

    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)
    #    fig, axs = plt.subplots(2, 2, hspace=1.0, wspace=1.0)
    (ax1, ax2), (ax3, ax4) = gs.subplots(sharex=False, sharey=False)
    ax1.plot(x, sypd, 'tab:blue')
    ax1.set(xlabel='leg', ylabel='SYPD')
    ax2.plot(x, elapsedtime, 'tab:orange')
    ax2.set(xlabel='leg', ylabel='Elapsed Time (sec)')
    ax3.plot(x, chpsy, 'tab:green')
    ax3.set(xlabel='leg', ylabel='CHPSY')
    ax4.plot(x, sbu, 'tab:red')
    ax4.set(xlabel='leg', ylabel='SBU')


def plot_sypd_vs_nptot(expnames, leg):

    nptot,sypd = group_by_nptot_sypd(expnames, leg)
    plt.xlabel("NPTOT")
    plt.ylabel("SYPD")
    pp = plt.scatter(nptot,sypd)
    pp = plt.plot(nptot,sypd)
    
    return pp


# multiplots
#  missing: averages and plots vs nptot e nfrac (for different expnames)
#  create plot at fixed number of nodes, ordered by nfrac
###############################################################

# Save table
# expnames = ['LB01', 'LB02', 'LB03', 'LB04', 'LB05', 'LB06', 'LB07', 'LB10', 'LB11', 'LB12', 'LB13', 'LB14', 'LB15', 'LB16', 'LB21', 'LB22', 'LB23', 'LB24', 'LB25', 'LB26', 'LB31', 'LB32', 'LB33', 'LB34', 'LB35', 'LB36', 'LB37', 'LC32' ]
def save_table(expnames, leg):

    npa=np.array([]); npo=np.array([])
    sypd=np.array([]); syph=np.array([])
    chpsy=np.array([]); elapsedtime=np.array([])
    sbu=np.array([])
    for expname in expnames:
        legs = selectleg(expname, leg)
        np1,np2 = compute_nodes(expname)
        npa = np.append(npa,np1)
        npo = np.append(npo,np2)
        sypd = np.append(sypd,np.mean(compute_sypd(expname, legs)))
        syph = np.append(syph,np.mean(compute_syph(expname, legs)))
        chpsy = np.append(chpsy,np.mean(compute_chpsy(expname, legs)))
        elapsedtime = np.append(elapsedtime,np.mean(compute_elapsedtime(expname, legs)))
        sbu = np.append(sbu,np.mean(compute_sbu(expname, legs)))

    nptot=npa+npo+1
    nfrac = npo/npa
    nodes = np.ceil(nptot/128).astype(int)    
    unique_nodes = set(nodes)
    sorted_nodes = sorted(unique_nodes)

    # sort vectors based on oifs/nemo fraction
    vars = [ nodes, nfrac, npa, npo, sypd, syph, chpsy, elapsedtime, sbu ]
    sorted_idx = np.argsort(nfrac)
    sorted_vars = [var[sorted_idx] for var in vars]
    sorted_expnames = [expnames[i] for i in sorted_idx]

    # write output
    with open('table_lobs.txt', 'w') as file:
        print('# Nodes(1) NPO/NPA(2) NPA(3) NPO(4) SYPD(5) SYPH(6) CHPSY(7) Elapsed_Time(8) SBU(9) expname(10) ', file=file)
        print('# ', file=file)
        for node in sorted_nodes:
            print('', file=file)
            print(f'# N={node}', file=file)
            for i in range(len(expnames)):
                if sorted_vars[0][i] == node: # checkes nodes[i]
                    row = [f"{var[i]:<5}" for var in sorted_vars]
                    row.append(f"({sorted_expnames[i]})")
                    print(" ".join(row), file=file)
        print('', file=file)
