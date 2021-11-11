import os, sys
import pathlib
from pathlib import Path
import datetime
import json
import re

import click
import numpy as np
import tifffile
import xmltodict
import yaml
import shutil

import dask
import dask.array as da

import pandas as pd

@click.group()
def cli():
    pass


@cli.command(name='get_xml_files')
@click.argument('filename', type=click.Path(exists=True, path_type=pathlib.Path))
def get_xml_files(folder):
    xml_files = sorted([file for file in folder.glob('**/Experiment.xml')
                        if 'zstack' not in file.parent.stem])
    return xml_files


def get_raw_files(folder):
    raw_files = sorted([file for file in folder.glob('**/Image_001_001.raw')
                        if 'zstack' not in file.parent.stem])
    return raw_files



def get_recursively(search_dict, field):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():

        if key == field:
            fields_found.append(value)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)

    return fields_found


def guess_type(s):
    if s == "":
        return None
    elif re.match("\A[0-9]+\.[0-9]+\Z", s):
        return float
    elif re.match("\A[0-9]+\Z", s):
        return int
    # 2019-01-01 or 01/01/2019 or 01/01/19
    elif re.match("\A[0-9]{4}-[0-9]{2}-[0-9]{2}\Z", s) or \
            re.match("\A[0-9]{2}/[0-9]{2}/([0-9]{2}|[0-9]{4})\Z", s):
        return datetime.date
    elif re.match("\A(true|false)\Z", s):
        return bool
    else:
        return str


def postprocessor(path, key, value):
    try:
        fcn = guess_type(value)
        return key, fcn(value)
    except (ValueError, TypeError):
        return key, value


# %%
# get_thorsync_xmlpaths
# self.sync_paths = sorted(glob.glob(self.sync_path))
# self.sync_episodes = [sorted(glob.glob(sync_path + "/Episode*.h5"))
#                       for sync_path in self.sync_paths]
# self.sync_xml = [sync_path + "/ThorRealTimeDataSettings.xml"
#                  for sync_path in self.sync_paths]

def read_metadata(filename):
    with open(filename) as f:
        ometa = xmltodict.parse(f.read(), postprocessor=postprocessor)
        meta = json.loads(json.dumps(ometa))
    return meta


def read_raw(filename):
    Yraw = np.fromfile(filename, dtype=np.uint16, count=- 1, offset=0)
    return Yraw


def get_utime(d):
    item = d['ThorImageExperiment']['Date']
    if isinstance(item, list):
        u_time = int(item[0]['@uTime'])
    elif isinstance(item, dict):
        u_time = int(item['@uTime'])
    return u_time


def get_px_size(d):
    lsm = d['ThorImageExperiment']['LSM']
    if isinstance(lsm, list):
        lsm = lsm[0]
    xpx = int(lsm['@pixelX'])
    ypx = int(lsm['@pixelX'])
    return xpx, ypx


def get_total_frames(d):
    n_frames = get_recursively(d, '@frames')[0]
    return n_frames


def get_timepoints(d):
    n_timepoints = get_recursively(d, '@timepoints')[0]
    return n_timepoints


def is_fast_z(d):
    fast_z = get_total_frames(d) == get_timepoints(d)
    return not fast_z


def is_zstack(d):
    T = get_timepoints(d)
    if T == 1:
        return True
    else:
        return False


def z_steps(d):
    if is_fast_z(d):
        steps = d['ThorImageExperiment']['ZStage'][0]['@steps']
    else:
        steps = 1
    return steps


def flyback_frames(d):
    if is_fast_z(d):
        ff = d['ThorImageExperiment']['Streaming'][0]['@flybackFrames']
    else:
        ff = 0
    return ff


def movie_shape(d):
    xpx, ypx = get_px_size(d)
    zpx = z_steps(d) + flyback_frames(d)
    T = get_timepoints(d)
    return T, zpx, ypx, xpx


def get_dimensions(xml_root):
    for child in xml_root:
        if child.tag == "LSM":
            xpx = int(child.attrib['pixelX'])
            ypx = int(child.attrib['pixelY'])
            if int(child.attrib['averageMode']) == 1:
                naverage = int(child.attrib['averageNum'])
            else:
                naverage = None
            if 'widthUM' in child.attrib:
                width_um = float(child.attrib['widthUM'])
                height_um = float(child.attrib['heightUM'])
    return xpx, ypx, width_um, height_um


def write_split_hyperstack(savedir, Y, sections=3, axis=0):
    for i, y in enumerate(np.split(Y, sections, axis=0)):
        fname = f"movie_{str(i).zfill(4)}.ome.tif"
        tifffile.imwrite(savedir.joinpath(fname), np.moveaxis(y.T, -1, 0), metadata=dict(axes='TZYX'))


def to_meta_dict(filename):
    meta = read_metadata(filename)

    meta_dict = dict(siz=list(movie_shape(meta)),
                     frames=get_total_frames(meta),
                     timepoints=get_timepoints(meta),
                     fast_z=is_fast_z(meta),
                     z_steps=z_steps(meta),
                     flyback_frames=flyback_frames(meta))
    return meta_dict


# @cli.command(name='xml_to_yaml')
# @click.argument('filename', type=click.Path(exists=True, path_type=pathlib.Path))
def xml_to_yaml(filename):
    meta_dict = to_meta_dict(filename)
    savename = filename.parent.joinpath('xml_meta.yaml')
    with open(savename, 'w') as f:
        yaml.dump(meta_dict, f)
    return savename


def xml_to_json(filename):
    meta_dict = to_meta_dict(filename)
    savename = filename.parent.joinpath('xml_meta.json')
    with open(savename, 'w') as f:
        json.dump(meta_dict, f, indent=4)
    return savename


@cli.command(name='test_click')
def test_click():
    click.echo('IT WORKED HOLY SHIT')


def dask_raw2tiffs(raw_file, n_blocks=None, drop_flybacks_frames=True):
    # load xml info for ThorImage movie
    yaml_file = raw_file.parent.joinpath("xml_meta.yaml")
    with open(yaml_file, 'r') as f:
        meta = yaml.safe_load(f)

    json_file = raw_file.parent.joinpath("sync_meta.json")
    with open(json_file, 'r') as f:
        sync_meta = json.load(f)

    # find # of blocks from xml_meta.yaml file
    if n_blocks is None:
        n_blocks = sync_meta['n_blocks']

    # get dimensions / chunking info
    T, *dims = meta['siz']
    d3, d2, d1 = dims
    chunks = int(np.prod(meta['siz'])/n_blocks)
    T_block = int(T/n_blocks)

    # files to save
    subdir = raw_file.parent.joinpath("raw_tifs")
    if not subdir.is_dir():
        subdir.mkdir()
    savenames = [subdir.joinpath(f"stk_{str(i).zfill(4)}.tif") for i in range(n_blocks)]

    # load raw image data
    print('reading raw file...')
    yraw = np.fromfile(raw_file, dtype='<u2')
    print('reshaping raw file...')
    yr = yraw.reshape(n_blocks, -1, d3, d2, d1)
    for i in range(n_blocks):
        y = yr[i, :, :meta['z_steps'], :, :]
        tifffile.imsave(savenames[i], y)
        print(f"{savenames[i]} saved successfully.")
        print(y.shape)
    return savenames


def move_to_processed_data(xml_meta_file):
    rel = xml_meta_file.relative_to(pathlib.Path.cwd(), 'data', '2p')
    to_file = pathlib.Path.cwd().joinpath('data', 'processed_data', rel)

    if not to_file.parent.is_dir():
        to_file.parent.mkdir(parents=True)

    shutil.copy(xml_meta_file, to_file)
    return to_file


def move_to_2p_data(xml_meta_file):
    rel = xml_meta_file.relative_to(pathlib.Path.cwd(), 'data', 'processed_data')
    to_file = pathlib.Path.cwd().joinpath('data', '2p', rel)
    shutil.copy(xml_meta_file, to_file)
    return to_file


def up_one_dir(filename):
    new_file = shutil.move(str(filename), str(filename.parents[1]))
    return new_file


def frametimeinfo_to_timestamps_npy(filename):
    df_fti = pd.read_csv(filename)
    if 'plane' in df_fti.columns:
        timestamps = df_fti.loc[df_fti['plane'] == 0, :]['frame_times'].to_numpy()
    else:
        timestamps = df_fti['frame_times'].to_numpy()
    savename = filename.parent.joinpath('timestamps.npy')
    np.save(savename, timestamps)

    return savename


if __name__ == '__main__':
    DATA_DIR = Path('./data/processed_data')

    file_list = sorted(list(DATA_DIR.rglob("2021-06-16/**/frametimeinfo.csv")))

    for i, f in enumerate(file_list):
        print(f"{i} ----- {f}")

    file_idx = int(input("Enter file #: (enter -1 to run all)"))
    if file_idx == -1:
        print("RUNNING ALL: ")
        print("------------")

    for f in file_list:
        print(frametimeinfo_to_timestamps_npy(f))

# DATA_2P_DIR = Path(
#     "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/2p/2021-08-24")
# DATA_2P_DIR = DATA_2P_DIR.relative_to(Path.cwd())
#
# thorsync_xml_files = sorted(list(DATA_2P_DIR.glob('**/Episode*.h5')))
#
# xml_filelist = sorted([file for file in DATA_2P_DIR.glob('**/Experiment.xml')
#                        if 'zstack' not in file.parent.stem])
#
# for xml_file in xml_filelist:
#     xml_to_yaml(xml_file)

# def save_mmap(savename, Y):
#     fp = np.memmap(savedir.joinpath("movie.mmap"), dtype='uint16', shape=Y0.T.shape, mode='w+')
#     fp[:] = Y0.T[:]
#     fp.flush()
#     newfp = np.memmap(savename, dtype='uint16', mode='r', shape=Y.T.shape)