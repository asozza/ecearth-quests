#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Command line tool to branch an EC-Earth4 experiment: create a new experiment
(with a new name) in the run directory, restarting from a given leg of an
existing experiment.

It is the "run-dir" counterpart of duplicate-job.py (which only duplicates the
job/launch directory). Typical workflow:

    python duplicate-job.py -i aa00 -o bb00          # job dir
    python branch_ece4.py aa00 bb00 12 --dry-run     # check what will happen
    python branch_ece4.py aa00 bb00 12               # run dir, restart from leg 12

What it does:
  1. copies the run directory of SRC into DST (symlinks to input data are kept
     as symlinks), skipping output/, restart/, log/, post/, saveic/, logs and
     the run-time restart files of the current leg;
  2. copies restart/LEG of SRC into DST/restart/LEG, renaming NEMO restart
     files from SRC_* to DST_* (rebuilt restart*.nc files are skipped);
  3. puts the restart files of LEG in place in the DST run dir (same logic as
     rollback_ece4: rstas/rstos/rcf copied, NEMO and srf files linked,
     srf renamed according to CTIME in rcf, time.step updated);
  4. updates leginfo.yml (experiment name, leg number, leg start/end date);
  5. replaces SRC with DST in the small text files of the run dir;
  6. optionally (--copy-output) copies NEMO/OIFS output preceding the branch
     point, so that the new experiment has a continuous time series.

The source experiment is never modified.

Authors: Alessandro Sozza (CNR-ISAC)
Date: Oct 2026
"""

import argparse
import fnmatch
import glob
import os
import re
import shutil

import yaml
from dateutil.relativedelta import relativedelta


# global paths
base_path = "/ec/res4/scratch/itas/ece4"

# files in the run dir that belong to the current leg and must not be copied
RUNTIME_FILES = ['rstas.nc', 'rstos.nc', 'srf*', 'restart*.nc', 'rcf',
                 '*_restart*.nc', 'time.step']
# folders of the run dir handled separately (or not copied at all)
SKIP_FOLDERS = ['output', 'restart', 'log', 'post', 'saveic']
# logs and leftovers
SKIP_MISC = ['*.log', 'ecsbatch*', 'ocean.output*', 'core*', '*~',
             'run.stat', 'layout.dat', 'output.namelist.*']

MAX_TEXT_SIZE = 1024 * 1024   # 1 MB: larger files are not edited

DRY_RUN = False


# ---------------------------------------------------------------------------
# small helpers that honour --dry-run
# ---------------------------------------------------------------------------

def log(msg):
    print(("[dry-run] " if DRY_RUN else "") + msg)


def mkdir(path):
    if not os.path.isdir(path):
        log(f"Creating folder {path}")
        if not DRY_RUN:
            os.makedirs(path)


def transfer(src, dst, hardlink=False):
    """Copy a file, or hardlink it (falls back to copy across filesystems)"""
    if DRY_RUN:
        log(f"{'Hardlinking' if hardlink else 'Copying'} {src} -> {dst}")
        return
    if hardlink:
        try:
            os.link(src, dst)
            return
        except OSError:
            pass
    shutil.copy2(src, dst)


def symlink(target, linkname):
    log(f"Linking {os.path.basename(linkname)} -> {target}")
    if not DRY_RUN:
        if os.path.lexists(linkname):
            os.remove(linkname)
        os.symlink(target, linkname)


def write_text(path, text):
    if not DRY_RUN:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(text)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def get_nemo_timestep(filename):
    """ Get timestep from a NEMO restart file (EXP_TIMESTEP_restart_XXXX.nc) """
    return os.path.basename(filename).split('_')[1]


def get_ctime_from_rcf(rcf_path):
    """ Extract CTIME value from the OIFS rcf namelist """
    if not os.path.exists(rcf_path):
        return None
    with open(rcf_path, 'r') as f:
        match = re.search(r'CTIME\s*=\s*"\s*(\d+)\s*"', f.read())
    return match.group(1).strip() if match else None


def rename_exp(name, src, dst):
    """ Replace the experiment name in a file name """
    return name.replace(src, dst)


def years_in_name(name):
    """ Plausible years found in date tokens (YYYY, YYYYMM, YYYYMMDD) """
    tokens = re.findall(r'(?<!\d)(\d{4})(?:\d{2}){0,2}(?!\d)', name)
    return [int(t) for t in tokens if 1000 <= int(t) <= 3000]


# ---------------------------------------------------------------------------
# main steps
# ---------------------------------------------------------------------------

def read_leginfo(path, src, dst, leg):
    """ Load leginfo.yml, rename the experiment and move it to the given leg """

    with open(path, 'r', encoding='utf-8') as f:
        text = f.read().replace(src, dst)
    leginfo = yaml.load(text, Loader=yaml.FullLoader)

    info = leginfo['base.context']['experiment']['schedule']['leg']
    if leg > info['num']:
        raise ValueError(f"Leg {leg} not reached yet by {src} "
                         f"(current leg is {info['num']}).")

    # leg length: from start/end if available, otherwise 1 year (as rollback)
    if 'end' in info and info['end'] is not None:
        step = relativedelta(info['end'], info['start'])
    else:
        step = relativedelta(years=1)

    delta = leg - info['num']
    newstart = info['start'] + step * delta
    log(f"leginfo.yml: leg {info['num']} -> {leg}, "
        f"start {info['start']} -> {newstart}")

    info['num'] = leg
    info['start'] = newstart
    if 'end' in info and info['end'] is not None:
        info['end'] = newstart + step

    return leginfo, newstart


def copy_rundir(src_dir, dst_dir, src):
    """ Copy the run directory, skipping outputs, restarts and logs """

    patterns = RUNTIME_FILES + SKIP_MISC + [f"{src}_*", f"ICM*{src}*"]

    def ignore(dirpath, names):
        if os.path.abspath(dirpath) != os.path.abspath(src_dir):
            return set()
        skipped = {n for n in names if n in SKIP_FOLDERS}
        for pat in patterns:
            skipped |= set(fnmatch.filter(names, pat))
        return skipped

    if DRY_RUN:
        names = sorted(os.listdir(src_dir))
        for n in sorted(set(names) - ignore(src_dir, names)):
            log(f"Copying run dir entry {n}")
        return

    log(f"Copying run dir {src_dir} -> {dst_dir}")
    shutil.copytree(src_dir, dst_dir, symlinks=True, ignore=ignore)

    # warn about absolute symlinks still pointing into the source experiment
    for root, dirs, files in os.walk(dst_dir):
        for n in dirs + files:
            p = os.path.join(root, n)
            if os.path.islink(p):
                target = os.readlink(p)
                if os.path.isabs(target) and target.startswith(src_dir + os.sep):
                    print(f"WARNING: {p} still points to {target}")


def copy_restart(src_rst, dst_rst, src, dst, hardlink):
    """ Copy restart/LEG, renaming NEMO restarts to the new experiment """

    mkdir(dst_rst)
    for file in sorted(glob.glob(os.path.join(src_rst, '*'))):
        base = os.path.basename(file)
        if not os.path.isfile(file):
            continue
        # skip rebuilt NEMO restarts (removed by rollback as well)
        if base.startswith('restart') and base.endswith('.nc'):
            continue
        newbase = dst + base[len(src):] if base.startswith(f"{src}_") else base
        transfer(file, os.path.join(dst_rst, newbase), hardlink)


def setup_restart(dst_dir, dst_rst, src_rst, src, dst):
    """ Put restart files of the branch leg in the run dir (as rollback_ece4) """

    # in dry-run the target restart folder does not exist: list the source one
    listing_dir = src_rst if DRY_RUN else dst_rst

    def target_name(base):
        return dst + base[len(src):] if (DRY_RUN and base.startswith(f"{src}_")) else base

    # OASIS restarts and rcf are copied
    for name in ['rstas.nc', 'rstos.nc', 'rcf']:
        if os.path.isfile(os.path.join(listing_dir, name)):
            log(f"Copying restart {name}")
            if not DRY_RUN:
                shutil.copy2(os.path.join(dst_rst, name), os.path.join(dst_dir, name))

    # NEMO restarts are linked: EXP_TIMESTEP_restart_XXXX.nc -> restart_XXXX.nc
    nemo_files = sorted(glob.glob(os.path.join(listing_dir, "*_restart*.nc")))
    for file in nemo_files:
        base = target_name(os.path.basename(file))
        linkname = os.path.join(dst_dir, '_'.join(base.split('_')[2:]))
        symlink(os.path.join(dst_rst, base), linkname)

    # time.step from NEMO restart
    if nemo_files:
        timestep = int(get_nemo_timestep(nemo_files[0]))
        log(f"Writing time.step = {timestep}")
        write_text(os.path.join(dst_dir, 'time.step'), str(timestep))
    else:
        print("WARNING: no NEMO restart found, time.step not written")

    # OIFS srf restarts are linked, renamed after CTIME in rcf
    ctime = get_ctime_from_rcf(os.path.join(src_rst, 'rcf'))
    for file in sorted(glob.glob(os.path.join(listing_dir, 'srf*'))):
        base = os.path.basename(file)
        ext = base.split('.')[-1]
        target = f"srf{ctime}.{ext}" if ctime else base
        symlink(os.path.join(dst_rst, base), os.path.join(dst_dir, target))


def rename_in_textfiles(dst_dir, src, dst, exclude=('leginfo.yml',)):
    """ Replace SRC with DST in small text files at the top of the run dir """

    for name in sorted(os.listdir(dst_dir)):
        path = os.path.join(dst_dir, name)
        if name in exclude or os.path.islink(path) or not os.path.isfile(path):
            continue
        if os.path.getsize(path) > MAX_TEXT_SIZE:
            continue
        try:
            with open(path, 'r', encoding='utf-8') as f:
                text = f.read()
        except (UnicodeDecodeError, OSError):
            continue
        if src in text:
            log(f"Replacing {src} -> {dst} in {name}")
            write_text(path, text.replace(src, dst))


def copy_output(src_dir, dst_dir, src, dst, branch_year, hardlink):
    """ Copy NEMO/OIFS output strictly preceding the branch year """

    src_out = os.path.join(src_dir, 'output')
    if not os.path.isdir(src_out):
        return
    ncopied = 0
    for root, _, files in os.walk(src_out):
        rel = os.path.relpath(root, src_out)
        outdir = os.path.join(dst_dir, 'output', rel)
        for name in sorted(files):
            years = years_in_name(name)
            if years and max(years) >= branch_year:
                continue
            if not DRY_RUN:
                os.makedirs(outdir, exist_ok=True)
            transfer(os.path.join(root, name),
                     os.path.join(outdir, rename_exp(name, src, dst)), hardlink)
            ncopied += 1
    print(f"Output files copied: {ncopied} (before {branch_year})")


def branch_ece4(src, dst, leg, copy_out=False, hardlink=False):
    """ Create experiment DST in the run dir, branching SRC at leg LEG """

    src_dir = os.path.join(base_path, src)
    dst_dir = os.path.join(base_path, dst)
    leg_str = str(leg).zfill(3)
    src_rst = os.path.join(src_dir, 'restart', leg_str)
    dst_rst = os.path.join(dst_dir, 'restart', leg_str)

    # sanity checks
    if not os.path.isdir(src_dir):
        raise ValueError(f"Source experiment {src_dir} does not exist.")
    if not os.path.isdir(src_rst):
        raise ValueError(f"Restart folder {src_rst} does not exist.")
    if os.path.exists(dst_dir):
        raise ValueError(f"Target experiment {dst_dir} already exists.")
    if not os.path.isfile(os.path.join(src_dir, 'leginfo.yml')):
        raise ValueError(f"leginfo.yml not found in {src_dir}.")

    # read and update leginfo first: fails early if leg is not valid
    leginfo, newstart = read_leginfo(os.path.join(src_dir, 'leginfo.yml'),
                                     src, dst, leg)

    print(f"Branching {src} -> {dst} at leg {leg} ({newstart})")

    # 1. run dir
    copy_rundir(src_dir, dst_dir, src)
    for sub in ['output/nemo', 'output/oifs', 'restart', 'log']:
        mkdir(os.path.join(dst_dir, sub))

    # 2. restart folder of the branch leg
    copy_restart(src_rst, dst_rst, src, dst, hardlink)

    # 3. restart files in the run dir
    setup_restart(dst_dir, dst_rst, src_rst, src, dst)

    # 4. leginfo.yml
    log("Writing leginfo.yml")
    if not DRY_RUN:
        with open(os.path.join(dst_dir, 'leginfo.yml'), 'w', encoding='utf8') as f:
            yaml.dump(leginfo, f, default_flow_style=False)

    # 5. experiment name in text files
    if not DRY_RUN:
        rename_in_textfiles(dst_dir, src, dst)

    # 6. output
    if copy_out:
        copy_output(src_dir, dst_dir, src, dst, newstart.year, hardlink)

    print(f"Done. Remember to duplicate the job folder as well "
          f"(duplicate-job.py -i {src} -o {dst}).")


def parse_args():
    """ Command line parser for branch_ece4 """

    parser = argparse.ArgumentParser(description="Branch an EC-Earth4 experiment from a given leg into a new experiment")
    parser.add_argument("src", metavar="SRC", help="Source experiment name (e.g. aa00)")
    parser.add_argument("dst", metavar="DST", help="New experiment name (e.g. bb00)")
    parser.add_argument("leg", metavar="LEG", type=int, help="Leg to restart from")
    parser.add_argument("--copy-output", action="store_true", help="Copy NEMO/OIFS output preceding the branch point")
    parser.add_argument("--hardlink", action="store_true",
                        help="Use hardlinks instead of copies for restart and output "
                             "files (saves space; files stay safe if the source is deleted)")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be done without touching anything")
    parser.add_argument("--base-path", default=base_path, help=f"Root of the run directories (default: {base_path})")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    if len(args.src) != 4 or len(args.dst) != 4:
        raise ValueError("Experiment names must be 4 characters long.")
    if args.src == args.dst:
        raise ValueError("Source and target experiment must be different.")

    DRY_RUN = args.dry_run
    base_path = args.base_path

    branch_ece4(args.src, args.dst, args.leg, copy_out=args.copy_output, hardlink=args.hardlink)