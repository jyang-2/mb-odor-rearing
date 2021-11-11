"""
Utilities for parsing fly #, dates, and other information from directory structure.

This was built for mb_odor_project, which is organized w/ <flydate>/<flynum>/<movie>/<block> structure.

"""
from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import namedtuple


AnalysisSet = namedtuple("AnalysisSet", ["flydate", "flynum", "mov", "is_block", "blk", "alg"])


def parse_dataset_info(filepath):
    """
    Given a filepath, returns metadata about the dataset in a namedtuple. (flydate, flynum, mov, is_blk, blk)

    Parameters
    ----------
    filepath : path to a file in the suite2p folder

    Returns
    -------
        AnalysisSet(flydate, flynum, mov, is_blk, blk): named tuple

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)

    flydate = parse_parent_date(filepath).name
    flynum = parse_parent_fly(filepath).name
    mov = parse_parent_movie(filepath).name
    is_blk = in_block_dir(filepath)
    if is_blk:
        blk = parse_parent_block(filepath).name
        blk_dir = parse_parent_block(filepath)
    else:
        blk = None
        blk_dir = None
    alg = analysis_type(filepath)

    alys_set = AnalysisSet(flydate, flynum, mov, is_blk, blk, alg)

    return alys_set


def is_date_dir(fname):
    """ Returns true if folder is a "date" directory (contains metadata.yaml file) """
    if isinstance(fname, str):
        fname = Path(fname)
    return fname.joinpath('metadata.yaml').is_file()


def is_fly_dir(fname):
    """
    Returns true if folder is a "fly"  (or "flynum") directory.

    It should contain flyinfo.xlsx, and the name should be a digit.
    """
    if isinstance(fname, str):
        fname = Path(fname)
    return fname.joinpath('flyinfo.xlsx').is_file() and str.isdigit(fname.name)


def is_movie_dir(fname):
    """ Returns true if folder is a "movie" directory ( must contain Experiment.xml) """
    if isinstance(fname, str):
        fname = Path(fname)
    return fname.joinpath('Experiment.xml').is_file()


def is_block_dir(fname):
    """
    Returns true if folder is a "block" directory.

    It should contain "metadata.csv", and the name should be a digit.

    Use ryim.chonk.split_metadata_xlsx to split it into metadata.csv files, if it can be split evenly.

    """
    if isinstance(fname, str):
        fname = Path(fname)
    return fname.joinpath('metadata.csv').is_file() and str.isdigit(fname.name)


def in_block_dir(fname):
    """Returns True if file is under a block directory"""
    if isinstance(fname, str):
        fname = Path(fname)
    for f in fname.parents:
        if is_block_dir(f):
            return True
    return False


def is_top_suite2p(fname):
    if isinstance(fname, str):
        fname = Path(fname)

    if fname.name == 'suite2p':
        return True
    else:
        return False


def parse_parent_suite2p(fname):
    if isinstance(fname, str):
        fname = Path(fname)

    for f in fname.parents:
        if is_top_suite2p(f):
            return f

    raise Exception("Not in a suite2p folder! ")


def parse_parent_block(fname):
    if isinstance(fname, str):
        fname = Path(fname)
    for f in fname.parents:
        if is_block_dir(f):
            return f
    return Path('')


def parse_parent_date(fname):
    if isinstance(fname, str):
        fname = Path(fname)
    for f in fname.parents:
        if is_date_dir(f):
            return f
    return None


def parse_parent_movie(fname):
    """ Get parent folder of any file, corresponding to entire movie folder """
    if isinstance(fname, str):
        fname = Path(fname)
    for f in fname.parents:
        if is_movie_dir(f):
            return f
    return None


# get parent folder of any file, corresponding to fly
def parse_parent_fly(fname):
    if isinstance(fname, str):
        fname = Path(fname)
    for f in fname.parents:
        if is_fly_dir(f):
            return f
    return None


def is_plane_folder(fname):
    """

    Parameters
    ----------
    fname : file path of directory.

    Returns
    -------
    True - if fname is a directory in a suite2p folder, and is named "combined', or plane<#>
    """
    if isinstance(fname, str):
        fname = Path(fname)

    if fname.is_dir() and ('suite2p' in fname.parent.name):
        if fname.name == 'combined':
            return True
        elif 'plane' in fname.name:
            return True
        else:
            return False
    else:
        return False


def analysis_type(fname):
    """
    Returns analysis type, {'suite2p', 'caiman', None}.
        - 'suite2p' if any of the parent folders contain 'suite2p' anywhere in the name
        - 'caiman' if any of the parent folders contain 'caiman' or anywhere in the name
        - None if neither condition is met


    Parameters
    ----------
    fname : filepath or directory

    Returns
    -------
        'suite2p' or 'caiman'

    """
    if isinstance(fname, str):
        fname = Path(fname)

    for f in fname.parents:
        folder = f.name
        caiman_flag = ("caiman" in folder) | ("cm_" in f.name) | ("_cm" in folder)
        suite2p_flag = ("suite2p" in folder) | ("s2p_" in folder) | ("_s2p" in folder)
        if suite2p_flag & ~caiman_flag:
            return "suite2p"
        elif ~suite2p_flag & caiman_flag:
            return "caiman"
        elif suite2p_flag & caiman_flag:
            return "CONFLICT"
    return None


def get_plane_number(fname):
    """

    Parameters
    ----------
    fname :

    Returns
    -------

    """
    if isinstance(fname, str):
        fname = Path(fname)

    if is_plane_folder(fname):
        if 'plane' in fname.name:
            return int(re.findall(r'\d+', fname.name)[0])
        elif fname.name == 'combined':
            subdirs = list(fname.parent.glob('plane*'))
            return sorted([get_plane_number(item) for item in subdirs])
    else:
        raise Exception("Not a plane folder! ")


def dir_info_str(filename):
    """
    Prints information abt dataset, based on path to file & directory structure.

    Parameters
    ----------
    filename : pathlib.Path, or string

    Returns
    -------

    """
    dir_info_str = f"""\nfilename = {filename},
        \n\tflydate: {parse_parent_date(filename).name}, 
        \n\tflynum: {parse_parent_fly(filename).name}",
        \n\tmovie: {parse_parent_movie(filename).name}",
        \n\tis in block folder: {in_block_dir(filename)},
        \n\tblock: {parse_parent_block(filename).name}"""

    return dir_info_str
